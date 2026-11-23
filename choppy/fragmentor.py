from dataclasses import dataclass, field, InitVar
import copy
from typing import NamedTuple
import warnings
import numpy as np
from choppy.constraints import Condition, MapStore, OverlapGrid, check_seq


class NoAssemblyError(RuntimeError):
    """Raised when no chain of overlaps spans the sequence."""


class GapError(NoAssemblyError):
    """Raised when stretches of the sequence lack candidate overlaps for longer than any
    segment can bridge.

    gaps holds the offending regions as (start, end), sorted by start, so callers can act
    on them (add candidates there, relax the constraints, or add local rules).
    """

    def __init__(self, message: str, gaps: list[tuple[int, int]]):
        super().__init__(message)
        self.gaps = gaps


class Candidates:
    """Candidate overlaps as parallel arrays, one entry per (start, end).

    Always sorted by start, then end, with no two entries sharing (start, end). The arrays
    are read-only so the order and the alignment of the fields cannot be broken from
    outside; every change goes through add, which re-sorts all fields together, or remove,
    which drops entries from all of them at once.

    The first and last entries are always the zero-length sentinels (0, 0) and
    (seq_len, seq_len) marking the sequence edges, so the stretches before the first real
    overlap and after the last one are ordinary gaps, and there is at least one gap even
    without any real overlaps. add cannot create or replace them: real overlaps need
    0 <= start < end <= seq_len.
    """

    def __init__(self, seq_len: int, starts=(), ends=(), extra_penalty=0.0):
        self.seq_len = seq_len
        self.starts = self._freeze(np.array([0, seq_len], dtype=np.int64))
        self.ends = self._freeze(np.array([0, seq_len], dtype=np.int64))
        self.extra_penalty = self._freeze(np.zeros(2, dtype=np.float64))
        self.add(starts, ends, extra_penalty)

    @staticmethod
    def _freeze(arr: np.ndarray) -> np.ndarray:
        arr.flags.writeable = False
        return arr

    def _pairs(self, starts, ends) -> tuple[np.ndarray, np.ndarray]:
        """starts and ends as 1-D int64 arrays of equal length, checked to be real overlaps:
        0 <= start < end <= seq_len, which the sentinels are not."""
        new_starts = np.atleast_1d(np.asarray(starts, dtype=np.int64))
        new_ends = np.atleast_1d(np.asarray(ends, dtype=np.int64))
        if new_starts.ndim != 1 or new_starts.shape != new_ends.shape:
            raise ValueError(
                f"starts and ends must be 1-D and of equal length, got shapes "
                f"{new_starts.shape} and {new_ends.shape}"
            )
        if np.any(new_starts < 0) or np.any(new_ends <= new_starts) or np.any(new_ends > self.seq_len):
            raise ValueError(f"Every candidate needs 0 <= start < end <= seq_len ({self.seq_len})")
        return new_starts, new_ends

    def _keys(self, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
        """One int per (start, end), equal exactly where both are."""
        return starts * (self.seq_len + 1) + ends

    def add(self, starts, ends, extra_penalty=0.0, skip_duplicates: bool = True) -> None:
        """Add candidates, keeping the order by start, then end.

        starts and ends are ints or array-likes of equal length; extra_penalty is a scalar
        for all of them or one value each. A candidate whose (start, end) is already
        present is ignored if skip_duplicates, otherwise its penalty replaces the stored
        one. Repeats within the added batch follow the same rule in the order given.
        """
        new_starts, new_ends = self._pairs(starts, ends)
        try:
            new_penalty = np.broadcast_to(np.asarray(extra_penalty, dtype=np.float64),
                                          new_starts.shape)
        except ValueError:
            raise ValueError(
                f"extra_penalty must be a scalar or have {len(new_starts)} values, "
                f"got shape {np.shape(extra_penalty)}"
            ) from None

        starts = np.concatenate((self.starts, new_starts))
        ends = np.concatenate((self.ends, new_ends))
        penalty = np.concatenate((self.extra_penalty, new_penalty))

        # lexsort is stable: among equal (start, end), stored entries come first, then the
        # new ones in the order given, so "first of a run" is the oldest, "last" the newest
        order = np.lexsort((ends, starts))
        starts, ends, penalty = starts[order], ends[order], penalty[order]

        repeat = (starts[1:] == starts[:-1]) & (ends[1:] == ends[:-1])
        keep = np.ones(len(starts), dtype=bool)
        if skip_duplicates:
            keep[1:] = ~repeat   # first of each run of equal (start, end)
        else:
            keep[:-1] = ~repeat  # last of each run

        self.starts = self._freeze(starts[keep])
        self.ends = self._freeze(ends[keep])
        self.extra_penalty = self._freeze(penalty[keep])

    def remove(self, starts, ends, missing_ok: bool = False) -> None:
        """Remove the candidates at (start, end), given as for add.

        One that is not present raises KeyError, and nothing is removed, unless missing_ok,
        which skips it. The sentinels cannot be removed: like add, remove takes only real
        overlaps, 0 <= start < end <= seq_len.
        """
        drop_starts, drop_ends = self._pairs(starts, ends)
        keys = self._keys(self.starts, self.ends)
        drop = self._keys(drop_starts, drop_ends)
        if not missing_ok:
            absent = ~np.isin(drop, keys)
            if absent.any():
                missing = list(zip(drop_starts[absent].tolist(), drop_ends[absent].tolist()))
                raise KeyError(f"No candidate at {missing}")
        keep = ~np.isin(keys, drop)
        self.starts = self._freeze(self.starts[keep])
        self.ends = self._freeze(self.ends[keep])
        self.extra_penalty = self._freeze(self.extra_penalty[keep])

    def __len__(self) -> int:
        """The number of entries, the two sentinels included."""
        return len(self.starts)

    @property
    def n_overlaps(self) -> int:
        """The number of real overlaps, without the sentinels."""
        return len(self) - 2

    def __iter__(self):
        """(start, end) pairs as plain ints, in order, the sentinels included."""
        return zip(self.starts.tolist(), self.ends.tolist())

    def __repr__(self) -> str:
        return f"Candidates({self.n_overlaps} overlaps)"


class _Generation(NamedTuple):
    """The last generation of a level's candidates (see Overlaps.generate_candidates)."""
    candidates: Candidates  # what it produced, with the overlaps added by hand merged in
    state: tuple            # what they depend on, as it was then (Overlaps._generation_state)
    own: Candidates         # what the generation itself produced, before that merge


@dataclass
class Overlaps:
    """Every potential overlap of a sequence, together with the conditions scored on it.

    The condition maps come from maps, a store that may be shared with other levels: its
    grid can be wider than this level's, and each condition reads a view of the shared
    map cut down to this level's grid. Without a store the level gets one of its own.

    Generated candidates remember the conditions and the min_step they were generated
    from, so the level can tell when they are out of date (see candidates_stale). Overlaps
    added by hand are kept apart in manual and merged back into every new generation,
    until remove_candidates takes them out again.

    min_step can be changed at any time; it is checked on assignment, and a step of 1,
    which thins nothing, is stored as None.
    """
    seq: str
    min_len: int # miminal overlap length
    max_len: int # maximal overlap length
    min_step: int | None = None # thinning parameter for overlaps (delta_start + delta_end >= min_step)

    conditions: list[Condition] = field(default_factory=list)
    candidates: Candidates | None = None  # None until generated
    maps: MapStore | None = field(default=None, repr=False)

    grid: OverlapGrid = field(init=False, repr=False)
    manual: Candidates = field(init=False, repr=False)  # the overlaps added by hand
    _generated: _Generation | None = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        self.grid = OverlapGrid(len(self.seq), self.min_len, self.max_len)
        self.manual = Candidates(len(self.seq))
        if self.maps is None:
            self.maps = MapStore(self.seq, self.grid)
        else:
            if self.maps.seq != self.seq:
                raise ValueError("The map store is for a different sequence")
            self.maps.grid.slices(self.grid)  # raises unless this level fits the store's grid

    def __setattr__(self, name: str, value) -> None:
        if name == "min_step":  # in __init__ too
            value = self._check_min_step(value)
        super().__setattr__(name, value)

    @staticmethod
    def _check_min_step(min_step) -> int | None:
        """min_step as stored: None, or an int of at least 2. A step of 1 thins nothing,
        so it is stored as None; a step below 1 is refused, as _thin cannot apply it."""
        if min_step is None:
            return None
        if isinstance(min_step, bool) or not isinstance(min_step, (int, np.integer)):
            raise TypeError(f"min_step must be a positive int or None, got {min_step!r}")
        if min_step <= 0:
            raise ValueError(f"min_step must be positive, got {min_step}")
        return None if min_step == 1 else int(min_step)

    @property
    def shape(self) -> tuple[int, int]:
        return self.grid.shape

    @property
    def constraints(self) -> list[Condition]:
        """The conditions with min_val and/or max_val, in the order added."""
        return [c for c in self.conditions if c.is_constraint]

    @property
    def objectives(self) -> list[Condition]:
        """The conditions with opt_val, in the order added."""
        return [c for c in self.conditions if c.is_objective]

    def add_condition(self, condition: Condition) -> None:
        """Bind condition to this level's view of its scorer's map and keep it.

        Its name, if set, must be unique here; that is checked before anything is bound, so
        a rejected condition is left untouched. Several conditions may share a scorer (and
        so its map), e.g. with different thresholds over different areas; all of them apply.
        """
        if condition.name is not None and any(c.name == condition.name for c in self.conditions):
            raise ValueError(f"There is already a condition named {condition.name!r} on this level")
        condition.bind(self.maps.view(condition.scorer, self.grid))
        self.conditions.append(condition)

    def get_condition(self, name: str) -> Condition:
        """The condition named name on this level."""
        for c in self.conditions:
            if c.name == name:
                return c
        names = [c.name for c in self.conditions if c.name is not None]
        raise KeyError(f"No condition named {name!r} on this level; named ones: {names}")

    def feasible_penalties(self) -> np.ndarray:
        """Extra penalty of each overlap under all constraints: the sum of their
        feasible() penalties, so 0 where every constraint passes outright, the relaxed
        penalties added up where some pass only relaxed, and inf where any forbids it."""
        penalties = np.zeros(self.shape, dtype=float)
        for constraint in self.constraints:
            penalties += constraint.feasible()
        return penalties

    def _thin(self, penalties: np.ndarray) -> np.ndarray:
        """Thin the feasible overlaps (finite penalties) to avoid excessive candidates.

        Keeps no two overlaps closer than min_step (|delta start| + |delta end|), setting
        the dropped ones to inf. The overlaps are taken greedily, smallest penalty first
        (then by start, then end), so a penalised overlap is kept only away from every
        kept better one. Works in place and returns penalties.
        """
        ## TODO (after objective scoring): thin with respect to objectives too
        if self.min_step is None:
            return penalties

        # all (i, -i + j) with |i| + |j| < min_step, before and after the overlap
        offsets = []
        for i in range(-self.min_step + 1, self.min_step):
            reach = self.min_step - abs(i)
            for j in range(-reach + 1, reach):
                if i == 0 and j == 0:
                    continue  # don't mask the overlap itself
                offsets.append((i, -i + j))
        offsets = np.array(offsets)

        for level in np.unique(penalties[~np.isinf(penalties)]):
            for i, j in np.argwhere(penalties == level):
                if np.isinf(penalties[i, j]):  # dropped by a kept overlap of this level
                    continue
                rows, cols = i + offsets[:, 0], j + offsets[:, 1]
                inside = (rows < penalties.shape[0]) & (cols < penalties.shape[1]) & (rows >= 0) & (cols >= 0)
                penalties[rows[inside], cols[inside]] = np.inf
        return penalties

    def generate_candidates(self, min_edge_len: int) -> Candidates:
        """Potential overlaps, sorted by start then end, between the two edge sentinels:
        the feasible ones, with the extra penalty of the relaxed constraints they pass, and
        every overlap added by hand, with its own (it replaces a generated one at the same
        start and end)."""
        penalties = self.feasible_penalties()
        penalties[~self.grid.edge_mask(min_edge_len)] = np.inf
        penalties = self._thin(penalties)

        rows, cols = np.nonzero(~np.isinf(penalties))
        starts, ends = self.grid.overlap(rows, cols)
        own = Candidates(len(self.seq), starts, ends, penalties[rows, cols])
        candidates = self._with_manual(own)

        self.candidates = candidates
        self._generated = _Generation(candidates, self._generation_state(), own)
        return candidates

    def _with_manual(self, own: Candidates) -> Candidates:
        """own with the overlaps added by hand merged in, as a separate Candidates; they
        replace any of own at the same start and end. own is left as it is."""
        merged = copy.copy(own)  # the arrays are read-only and add replaces them, so sharing is safe
        m = self.manual
        if m.n_overlaps:
            merged.add(m.starts[1:-1], m.ends[1:-1], m.extra_penalty[1:-1], skip_duplicates=False)
        return merged

    def _generation_state(self) -> tuple:
        """What the candidates depend on: min_step, and each condition with the version of
        its thresholds."""
        return self.min_step, tuple((c, c.version) for c in self.conditions)

    @property
    def candidates_stale(self) -> bool:
        """Whether the candidates have to be (re)generated before use.

        True when there are none yet, or when they were generated and what they depend on
        has changed since: a condition added, the thresholds of one updated, or min_step
        changed. Candidates assigned directly rather than generated are taken as given and
        never stale.
        """
        if self.candidates is None:
            return True
        if self._generated is None or self._generated.candidates is not self.candidates:
            return False
        return self._generated.state != self._generation_state()

    def add_candidates(self, starts, ends, extra_penalty=0.0) -> None:
        """Add overlaps by hand, e.g. to bridge a gap the constraints leave.

        They are kept in manual and survive every regeneration of the candidates (until
        remove_candidates), and are added to the current candidates too, if there are any.
        A hand-added overlap replaces whatever is at the same start and end, an earlier
        hand-added one included, so its extra_penalty is the one that applies.
        """
        self.manual.add(starts, ends, extra_penalty, skip_duplicates=False)
        if self.candidates is not None:
            self.candidates.add(starts, ends, extra_penalty, skip_duplicates=False)

    def remove_candidates(self, starts, ends) -> None:
        """Remove overlaps added by hand (see add_candidates), given as for add_candidates.

        They are dropped from manual, so later generations no longer bring them back, and
        from the current candidates. Where one had replaced a generated candidate at the
        same start and end, that one is back, with its own extra_penalty. Removing does not
        make the candidates stale. Raises KeyError, and changes nothing, if one of them was
        not added by hand.
        """
        self.manual.remove(starts, ends)
        if self.candidates is None:
            return
        generated = self._generated
        if generated is not None and generated.candidates is self.candidates:
            self.candidates = self._with_manual(generated.own)
            self._generated = generated._replace(candidates=self.candidates)
        else:  # assigned directly: taken as given, so only the removed overlaps go
            self.candidates.remove(starts, ends, missing_ok=True)

    def _require_candidates(self) -> Candidates:
        if self.candidates is None:
            raise ValueError("No candidates have been generated yet")
        return self.candidates

    def _gaps(self) -> tuple[np.ndarray, np.ndarray]:
        """The shortest segment across each split between consecutive candidates, as
        (opens, closes) arrays, one entry per split.

        A chain of overlaps runs forward through the candidates, so it crosses the split
        after candidate i with a segment opening on some candidate at or before i and
        closing on one after it. Starts are sorted, so the latest opening is start_i; ends
        are not, so the earliest closing is the minimum over all later ends, not end_{i+1}.
        The sentinels make the first and last splits run from and to the sequence edges.
        """
        c = self._require_candidates()
        closes = np.minimum.accumulate(c.ends[::-1])[::-1][1:]
        return c.starts[:-1], closes

    def gap_lengths(self) -> np.ndarray:
        """The length of the shortest segment across each split between candidates."""
        opens, closes = self._gaps()
        return closes - opens

    def plot_gap_lengths(self, min_gap_len: int = None, ax=None):
        """Histogram of the shortest segment lengths across the splits between candidates."""
        lengths = self.gap_lengths()
        if min_gap_len is not None:
            lengths = lengths[lengths >= min_gap_len]
        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
        ax.hist(lengths, bins=range(1, lengths.max(initial=0) + 2), align="left", color="tab:blue")
        ax.set_xlabel("gap length")
        ax.set_ylabel("count")
        ax.set_title(f"gaps between {self.candidates.n_overlaps} candidate overlaps")
        return ax

    def plot_rows(self, indices: int | str | list[int | str] | None = None, a: int = 0, b: int = None):
        """Plot the rows of several constraints as a column of subplots sharing the x axis.

        indices selects constraints by their position in self.constraints (all by default)
        or conditions by name, mixed freely. The overlap start range [a, b) is the same for
        every subplot. Returns the axes array.
        """
        if indices is None:
            indices = range(len(self.constraints))
        elif isinstance(indices, (int, str)):
            indices = [indices]
        constraints = self.constraints
        selected = []  # (subplot title, condition)
        for i in indices:
            if isinstance(i, str):
                c = self.get_condition(i)
                selected.append((c.label, c))
            else:
                c = constraints[i]
                selected.append((f"[{i}] {c.label}", c))
        if not selected:
            raise ValueError("There are no constraints to plot")
        if b is None:
            b = self.shape[0]

        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(len(selected), 1, sharex=True, squeeze=False,
                                 figsize=(8, 2.5 * len(selected)))
        axes = axes[:, 0]
        for ax, (title, condition) in zip(axes, selected):
            condition.plot_rows(a, b, ax=ax)
            ax.set_title(title)
            ax.set_xlabel("")
        axes[-1].set_xlabel("overlap start")
        fig.suptitle(f"Overlaps starting in [{a}, {b})")
        fig.tight_layout()
        return axes

    def get_gapped_regions(self, min_gap_len: int) -> list[tuple[int, int]]:
        """The shortest segments across the splits between candidates that are at least
        min_gap_len long.

        Each region runs from the start of the overlap opening the segment to the end of
        the one closing it, both overlaps included, so its length compares directly with
        the segment lengths. Returns a list of (start, end) tuples, sorted by start; the
        regions of neighbouring splits may overlap.
        """
        opens, closes = self._gaps()
        keep = closes - opens >= min_gap_len
        return list(zip(opens[keep].tolist(), closes[keep].tolist()))


@dataclass
class Assembly:
    """One cut of the sequence: the chosen overlaps and the segments between them."""
    segments: list[tuple[int, int]]
    overlaps: list[tuple[int, int]]

    def __len__(self) -> int:
        return len(self.segments)


@dataclass
class Fragmentor:
    seq: str
    min_segment_overlap_len: InitVar[int]
    max_segment_overlap_len: InitVar[int]
    min_segment_len: int
    max_segment_len: int
    min_fragment_overlap_len: InitVar[int | None] = None
    max_fragment_overlap_len: InitVar[int | None] = None
    min_fragment_len: int | None = None
    max_fragment_len: int | None = None
    # what both levels start with; each level's can be changed later on its own
    # (segment_overlaps.min_step, fragment_overlaps.min_step), and this one does not follow
    min_step: int | None = None

    grid: OverlapGrid = field(init=False, repr=False)  # spans the overlap lengths of both levels
    maps: MapStore = field(init=False, repr=False)     # shared by both levels
    fragment_overlaps: Overlaps | None = field(init=False, repr=False)
    segment_overlaps: Overlaps = field(init=False, repr=False)

    def __post_init__(self, min_segment_overlap_len, max_segment_overlap_len,
                 min_fragment_overlap_len, max_fragment_overlap_len):
        self.seq = check_seq(self.seq.upper())
        self.seq_len = len(self.seq)

        Overlaps._check_min_step(self.min_step)  # early, before the lengths; the levels check it again

        max_segment_overlap_len = self._check_lengths(
            "segment",
            min_segment_overlap_len,
            max_segment_overlap_len,
            self.min_segment_len,
            self.max_segment_len,
        )

        fragment_params = (
            min_fragment_overlap_len,
            max_fragment_overlap_len,
            self.min_fragment_len,
            self.max_fragment_len,
        )
        n_set = sum(p is not None for p in fragment_params)
        if n_set not in (0, 4):
            raise ValueError("Fragment lengths must be either all set or all None")
        self.use_fragments = n_set == 4
        if self.use_fragments:
            max_fragment_overlap_len = self._check_lengths("fragment", *fragment_params)
            if self.min_fragment_len > self.min_segment_len:
                raise ValueError(
                    f"min_fragment_len ({self.min_fragment_len}) must be <= min_segment_len "
                    f"({self.min_segment_len}): a fragment is a piece of a segment"
                )
            if self.max_fragment_len > self.max_segment_len:
                warnings.warn(
                    f"max_fragment_len ({self.max_fragment_len}) is longer than "
                    f"max_segment_len ({self.max_segment_len}): a fragment is a piece of a "
                    f"segment. Lowering it to {self.max_segment_len}.",
                    stacklevel=3,
                )
                # still > min_fragment_len, which is <= min_segment_len < max_segment_len
                self.max_fragment_len = self.max_segment_len

        if self.seq_len < self.min_segment_len:
            raise ValueError(
                f"Sequence is {self.seq_len} bp, shorter than min_segment_len "
                f"({self.min_segment_len}): there is nothing to cut"
            )

        ranges = [(min_segment_overlap_len, max_segment_overlap_len)]
        if self.use_fragments:
            ranges.append((min_fragment_overlap_len, max_fragment_overlap_len))
        self.grid = OverlapGrid(self.seq_len, min(lo for lo, _ in ranges), max(hi for _, hi in ranges))
        self._warn_unused_lengths(ranges)
        self.maps = MapStore(self.seq, self.grid)

        self.fragment_overlaps = Overlaps(self.seq, min_fragment_overlap_len, max_fragment_overlap_len,
                                          self.min_step, maps=self.maps) if self.use_fragments else None
        self.segment_overlaps = Overlaps(self.seq, min_segment_overlap_len, max_segment_overlap_len,
                                         self.min_step, maps=self.maps)

    def _warn_unused_lengths(self, ranges: list[tuple[int, int]]) -> None:
        """Warn when the shared grid is mostly overlap lengths no level uses.

        The two levels share one map per scorer over every length from the shortest
        minimum to the longest maximum. When their ranges are far apart, the lengths in
        between are scored and stored for nothing; past the point where they outnumber
        the used ones, two separate maps would have been smaller than the shared one.
        """
        used = len(set().union(*(range(lo, hi + 1) for lo, hi in ranges)))
        unused = self.grid.shape[1] - used
        if unused > used:
            warnings.warn(
                f"The segment and fragment overlap lengths {ranges} are far apart: the maps "
                f"shared by both levels span {self.grid.min_len}-{self.grid.max_len} bp, and "
                f"{unused} of those {self.grid.shape[1]} lengths are used by neither level",
                stacklevel=4,
            )

    @staticmethod
    def _check_lengths(name: str, min_overlap_len, max_overlap_len, min_len, max_len) -> int:
        """Validate one level's four lengths; return max_overlap_len, lowered if needed.

        min_len > 2 * max_overlap_len is the load-bearing one: besides being the geometric
        requirement, it is what lets OverlapGrid.edge_mask(min_len) stand in for a separate
        check that an overlap ends inside the sequence, so conditions never have to mask
        off the cells that run past the end.
        """
        values = {
            f"min_{name}_overlap_len": min_overlap_len,
            f"max_{name}_overlap_len": max_overlap_len,
            f"min_{name}_len": min_len,
            f"max_{name}_len": max_len,
        }
        for key, val in values.items():
            if val is None:
                raise ValueError(f"{key} must be set")
            if not isinstance(val, int) or isinstance(val, bool):
                raise ValueError(f"{key} must be int, got {type(val).__name__}")
            if val <= 0:
                raise ValueError(f"{key} must be positive, got {val}")
        if min_overlap_len > max_overlap_len:
            raise ValueError(
                f"min_{name}_overlap_len ({min_overlap_len}) > max_{name}_overlap_len ({max_overlap_len})"
            )
        if max_len <= min_len:
            raise ValueError(f"max_{name}_len ({max_len}) must be > min_{name}_len ({min_len})")
        if min_len <= 2 * max_overlap_len:
            longest = (min_len - 1) // 2  # the largest m with min_len > 2 * m
            if min_overlap_len > longest:
                raise ValueError(
                    f"min_{name}_overlap_len ({min_overlap_len}) leaves no room in a "
                    f"min_{name}_len ({min_len}) {name}: a {name} has to hold an overlap at "
                    f"each end without the two meeting, so no overlap may be longer than "
                    f"{longest}"
                )
            warnings.warn(
                f"max_{name}_overlap_len ({max_overlap_len}) is too long for "
                f"min_{name}_len ({min_len}): a {name} has to hold a maximal overlap at "
                f"each end without the two meeting. Lowering it to {longest}.",
                stacklevel=4,
            )
            return longest
        return max_overlap_len

    def _level(self, fragment: bool) -> Overlaps:
        """The fragment or the segment level."""
        if not fragment:
            return self.segment_overlaps
        if self.fragment_overlaps is None:
            raise ValueError("Fragment lengths were not set, so there is no fragment level")
        return self.fragment_overlaps

    def add_condition(self, condition: Condition, fragment: bool = False) -> None:
        self._level(fragment).add_condition(condition)

    def get_condition(self, name: str, fragment: bool = False) -> Condition:
        return self._level(fragment).get_condition(name)

    def _lengths(self, fragment: bool) -> tuple[int, int]:
        """min_len and max_len of the fragments or the segments."""
        self._level(fragment)  # raises if there is no fragment level
        if fragment:
            return self.min_fragment_len, self.max_fragment_len
        return self.min_segment_len, self.max_segment_len

    def generate_candidates(self, fragment: bool = False) -> Candidates:
        """Generate the level's candidates afresh, whether or not they are stale."""
        return self._level(fragment).generate_candidates(self._lengths(fragment)[0])

    def _current_candidates(self, fragment: bool) -> Candidates:
        """The level's candidates, regenerated first if they are stale."""
        level = self._level(fragment)
        if level.candidates_stale:
            self.generate_candidates(fragment)
        return level.candidates

    def plot_gap_lengths(self, fragment: bool = False, min_gap_len: int = None, ax=None):
        self._current_candidates(fragment)
        return self._level(fragment).plot_gap_lengths(min_gap_len, ax)

    def plot_rows(self, fragment: bool = False, indices: int | str | list[int | str] | None = None,
                  a: int = 0, b: int = None):
        return self._level(fragment).plot_rows(indices, a, b)

    def gap_check(self, fragment: bool = False, min_gap_len: int = None) -> None:
        """Raise GapError if some stretch between candidate overlaps is too long to bridge.

        A gap is the shortest segment that can cross a split between consecutive
        candidates: from the start of the overlap opening it to the end of the one closing
        it (see Overlaps.get_gapped_regions), with the sequence edges as zero-length
        overlaps. Every assembly crosses every split, so a gap longer than max_len rules
        out every assembly, whatever the rest of the chain looks like. (Fragments alike,
        with the fragment lengths.)

        Passing min_gap_len replaces that rule: every gap at least min_gap_len long is
        reported.

        Stale candidates (none yet, or conditions changed since) are regenerated first.
        """
        name = "fragment" if fragment else "segment"
        overlaps, max_len = self._level(fragment), self._lengths(fragment)[1]
        self._current_candidates(fragment)

        threshold = max_len + 1 if min_gap_len is None else min_gap_len
        gaps = overlaps.get_gapped_regions(threshold)

        if not gaps:
            return

        if min_gap_len is None:
            why = (f"longer than any {name} of at most {max_len} bp, so no assembly "
                   f"exists")
        else:
            why = f"at least {min_gap_len} bp long"
        lines = "\n".join(f"  [{a}, {b}): {b - a} bp" for a, b in gaps)
        raise GapError(
            f"{len(gaps)} region(s) of the sequence can only be crossed by a {name} "
            f"{why} (each region runs from the start of the candidate overlap before it "
            f"to the end of the one after it, or from/to the sequence edge):\n{lines}\n"
            f"To get past them, for each region either add candidate overlaps manually, "
            f"relax the existing constraints, or add relaxed local rules for the affected "
            f"region.",
            gaps,
        )

    def add_candidates(self, starts, ends, extra_penalty=0.0, fragment: bool = False) -> None:
        """Add overlaps by hand; they survive regeneration (see Overlaps.add_candidates)."""
        self._level(fragment).add_candidates(starts, ends, extra_penalty)

    def remove_candidates(self, starts, ends, fragment: bool = False) -> None:
        """Remove overlaps added by hand (see Overlaps.remove_candidates)."""
        self._level(fragment).remove_candidates(starts, ends)

    def explore(self, fragment: bool = False, show: int | str | list[int | str] | None = None,
                a: int = 0, b: int | None = None, row_height: int = 280):
        """An interactive editor of the level's conditions and candidates: a notebook widget,
        which needs the `gui` extra (see choppy.explorer).

        Each shown condition gets a plot of its row minima, as plot_rows draws them, with
        its thresholds, the bounds of its relaxed versions over the starts they apply to,
        and the candidates with the gaps between them. Next to the plot are its thresholds,
        area and relaxed versions, all editable. Edits act on this Fragmentor directly;
        the candidates are regenerated only by the button, never by an edit.

        show picks the conditions shown at first, by position in the level's conditions
        (objectives included, unlike plot_rows) or by name, mixed freely; the first four by
        default, and the rest can be ticked in the list. [a, b) is the range of overlap
        starts in view at first (all by default); every start is plotted, so the view can
        be zoomed out from there. row_height is the height in pixels of each condition's
        plot and panel.
        """
        level = self._level(fragment)  # raises if there is no fragment level, extra or not
        from choppy.explorer import LevelExplorer
        return LevelExplorer(level, self._lengths(fragment), name="fragment" if fragment else "segment",
                             show=show, a=a, b=b, row_height=row_height)

    def chop(self) -> Assembly:
        """Cut the sequence into segments at candidate overlaps, at the lowest cost.

        The cost is the number of segments plus the extra_penalty of every overlap used,
        so with no penalties it is the fewest segments, and a penalty of p is worth p
        extra segments: an overlap that passes a constraint only relaxed is used when
        avoiding it would cost more than that.

        Candidates that are missing or stale are (re)generated first, on both levels.
        """
        self.gap_check(fragment=False)
        if self.fragment_overlaps is not None:
            self.gap_check(fragment=True)

        c = self.segment_overlaps.candidates
        starts, ends = c.starts, c.ends

        weights, back = self._shortest_path(starts, ends, self.min_segment_len, self.max_segment_len,
                                            c.extra_penalty)
        if not np.isfinite(weights[-1]):
            raise NoAssemblyError(self._why_stuck(starts, ends, weights))

        path = self._walk_back(back)
        return Assembly(
            segments=[(int(starts[i]), int(ends[j])) for i, j in zip(path, path[1:])],
            overlaps=[(int(starts[i]), int(ends[i])) for i in path[1:-1]],
        )

    @staticmethod
    def _shortest_path(starts: np.ndarray, ends: np.ndarray, min_len: int, max_len: int,
                       extra_penalty: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
        """Cheapest path over the overlaps, each edge one segment.

        The overlaps are sorted by start and every edge runs forward, so one sweep in that
        order settles each vertex before it is used. An edge from i to j is the segment
        [start_i, end_j); it carries weight 1 plus the extra_penalty of overlap j (0 for
        the end sentinel), so a path costs its segments plus the penalties of the overlaps
        it uses, and without penalties a shortest path is a shortest assembly.

        Each vertex relaxes its whole window at once. Starts are sorted, so the window ends
        at the first start out of segment reach, and the segment length is then checked
        against both bounds per candidate.
        """
        n = len(starts)
        if extra_penalty is None:
            extra_penalty = np.zeros(n)
        weights = np.full(n, np.inf)
        back = np.full(n, -1, dtype=np.int64)
        weights[0] = 0.0

        for i in range(n):
            if not np.isfinite(weights[i]):
                continue
            opened = starts[i]
            stop = int(np.searchsorted(starts, opened + max_len, side="right"))
            if stop <= i + 1:
                continue

            window = slice(i + 1, stop)
            length = ends[window] - opened
            weight = weights[i] + 1.0 + extra_penalty[window]
            better = (length >= min_len) & (length <= max_len) & (weight < weights[window])

            improved = np.flatnonzero(better)
            weights[improved + i + 1] = weight[improved]
            back[improved + i + 1] = i

        return weights, back

    @staticmethod
    def _walk_back(back: np.ndarray) -> list[int]:
        path = [len(back) - 1]
        while back[path[-1]] != -1:
            path.append(int(back[path[-1]]))
        return path[::-1]

    def _why_stuck(self, starts: np.ndarray, ends: np.ndarray, weights: np.ndarray) -> str:
        """Report where the chain of overlaps ran out."""
        n_candidates = len(starts) - 2
        if n_candidates == 0:
            return ("No candidate overlaps at all: every potential overlap was rejected by the "
                    "constraints. Loosen the thresholds or widen the overlap length range.")

        reached = np.flatnonzero(np.isfinite(weights))
        frontier = int(ends[reached].max())
        after = starts[1:-1][starts[1:-1] > frontier]
        nearest = int(after.min()) if after.size else self.seq_len

        return (
            f"No assembly spans the sequence. {n_candidates} candidate overlaps, "
            f"the chain reaches position {frontier} of {self.seq_len} and stops there: "
            f"the next usable overlap starts at {nearest}, "
            f"{nearest - frontier} bp further on, which no segment of "
            f"{self.min_segment_len}-{self.max_segment_len} bp can bridge."
        )
