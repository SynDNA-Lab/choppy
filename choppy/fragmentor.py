from dataclasses import dataclass, field, InitVar
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
    outside; every change goes through add, which re-sorts all fields together.

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

    def add(self, starts, ends, extra_penalty=0.0, skip_duplicates: bool = True) -> None:
        """Add candidates, keeping the order by start, then end.

        starts and ends are ints or array-likes of equal length; extra_penalty is a scalar
        for all of them or one value each. A candidate whose (start, end) is already
        present is ignored if skip_duplicates, otherwise its penalty replaces the stored
        one. Repeats within the added batch follow the same rule in the order given.
        """
        new_starts = np.atleast_1d(np.asarray(starts, dtype=np.int64))
        new_ends = np.atleast_1d(np.asarray(ends, dtype=np.int64))
        if new_starts.ndim != 1 or new_starts.shape != new_ends.shape:
            raise ValueError(
                f"starts and ends must be 1-D and of equal length, got shapes "
                f"{new_starts.shape} and {new_ends.shape}"
            )
        if np.any(new_starts < 0) or np.any(new_ends <= new_starts) or np.any(new_ends > self.seq_len):
            raise ValueError(f"Every candidate needs 0 <= start < end <= seq_len ({self.seq_len})")
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


@dataclass
class Overlaps:
    """Every potential overlap of a sequence, together with the conditions scored on it.

    The condition maps come from maps, a store that may be shared with other levels: its
    grid can be wider than this level's, and each condition reads a view of the shared
    map cut down to this level's grid. Without a store the level gets one of its own.
    """
    seq: str
    min_len: int # miminal overlap length
    max_len: int # maximal overlap length
    min_step: int | None = None # thinning parameter for overlaps (delta_start + delta_end >= min_step)

    objectives: list[Condition] = field(default_factory=list)
    constraints: list[Condition] = field(default_factory=list)
    candidates: Candidates | None = None  # None until generated
    maps: MapStore | None = field(default=None, repr=False)

    grid: OverlapGrid = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self.grid = OverlapGrid(len(self.seq), self.min_len, self.max_len)
        if self.maps is None:
            self.maps = MapStore(self.seq, self.grid)
        else:
            if self.maps.seq != self.seq:
                raise ValueError("The map store is for a different sequence")
            self.maps.grid.slices(self.grid)  # raises unless this level fits the store's grid
        if self.min_step == 1:
            self.min_step = None  # no thinning needed if min_step is 1

    @property
    def shape(self) -> tuple[int, int]:
        return self.grid.shape

    def add_condition(self, condition: Condition) -> None:
        condition.bind(self.maps.view(condition.scorer, self.grid))
        if condition.opt_val is not None:
            self.objectives.append(condition)
        if condition.min_val is not None or condition.max_val is not None:
            self.constraints.append(condition)

    def feasible_mask(self) -> np.ndarray:
        """Boolean mask of the overlaps that satisfy all constraints."""

        mask = np.ones(self.shape, dtype=bool)
        for constraint in self.constraints:
            mask &= constraint.feasible()
        return mask

    def _thin(self, mask: np.ndarray) -> np.ndarray:
        ## TODO (after objective scoring): thin with respect to objectives, keeping the best options
        """Thins the mask of feasible overlaps to avoid excessive candidates."""

        if self.min_step is None:
            return mask

        offsets = []
        for i in range(self.min_step):
            for j in range(self.min_step - i):
                offsets.append((i, -i + j))
                if j > 0:
                    offsets.append((i, -i - j))
                    
        offsets.pop(0)  # don't mask the overlap itself
        offsets = np.array(offsets)

        for i, j in np.argwhere(mask):
            if not mask[i, j]:
                continue

            rows, cols = i + offsets[:, 0], j + offsets[:, 1]
            inside = (rows < mask.shape[0]) & (cols < mask.shape[1]) & (rows >= 0) & (cols >= 0)
            mask[rows[inside], cols[inside]] = False
        return mask

    def generate_candidates(self, min_edge_len: int) -> Candidates:
        """Potential overlaps, sorted by start then end, with no extra penalty, between the
        two edge sentinels."""
        mask = self.feasible_mask() & self.grid.edge_mask(min_edge_len)
        mask = self._thin(mask)

        starts, ends = self.grid.overlap(*np.nonzero(mask))
        self.candidates = Candidates(len(self.seq), starts, ends)

        return self.candidates

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

    def plot_rows(self, indices: int | list[int] | None = None, a: int = 0, b: int = None):
        """Plot the rows of several constraints as a column of subplots sharing the x axis.

        indices selects constraints by their position in self.constraints (all by default).
        The overlap start range [a, b) is the same for every subplot. Returns the axes array.
        """
        if indices is None:
            indices = range(len(self.constraints))
        elif isinstance(indices, int):
            indices = [indices]
        selected = [(i, self.constraints[i]) for i in indices]
        if not selected:
            raise ValueError("There are no constraints to plot")
        if b is None:
            b = self.shape[0]

        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(len(selected), 1, sharex=True, squeeze=False,
                                 figsize=(8, 2.5 * len(selected)))
        axes = axes[:, 0]
        for ax, (i, constraint) in zip(axes, selected):
            constraint.plot_rows(a, b, ax=ax)
            ax.set_title(f"[{i}] {type(constraint.scorer).__name__}")
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
    min_step: int | None = None

    grid: OverlapGrid = field(init=False, repr=False)  # spans the overlap lengths of both levels
    maps: MapStore = field(init=False, repr=False)     # shared by both levels
    fragment_overlaps: Overlaps | None = field(init=False, repr=False)
    segment_overlaps: Overlaps = field(init=False, repr=False)

    def __post_init__(self, min_segment_overlap_len, max_segment_overlap_len,
                 min_fragment_overlap_len, max_fragment_overlap_len):
        self.seq = check_seq(self.seq.upper())
        self.seq_len = len(self.seq)

        if self.min_step is not None and self.min_step <= 0:
            raise ValueError(f"min_step must be positive, got {self.min_step}")

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

    def add_condition(self, condition: Condition, fragment: bool = False) -> None:
        if fragment:
            if self.fragment_overlaps is None:
                raise ValueError("Fragment lengths were not set, so there is no fragment level")
            self.fragment_overlaps.add_condition(condition)
        else:
            self.segment_overlaps.add_condition(condition)

    def generate_candidates(self, fragment: bool = False) -> Candidates:
        if fragment:
            if self.fragment_overlaps is None:
                raise ValueError("Fragment lengths were not set, so there is no fragment level")
            return self.fragment_overlaps.generate_candidates(self.min_fragment_len)
        else:
            return self.segment_overlaps.generate_candidates(self.min_segment_len)

    def plot_gap_lengths(self, fragment: bool = False, min_gap_len: int = None, ax=None):
        if fragment:
            if self.fragment_overlaps is None:
                raise ValueError("Fragment lengths were not set, so there is no fragment level")
            return self.fragment_overlaps.plot_gap_lengths(min_gap_len, ax)
        else:
            return self.segment_overlaps.plot_gap_lengths(min_gap_len, ax)

    def plot_rows(self, fragment: bool = False, indices: int | list[int] | None = None,
                  a: int = 0, b: int = None):
        if fragment:
            if self.fragment_overlaps is None:
                raise ValueError("Fragment lengths were not set, so there is no fragment level")
            return self.fragment_overlaps.plot_rows(indices, a, b)
        else:
            return self.segment_overlaps.plot_rows(indices, a, b)

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
        """
        name = "fragment" if fragment else "segment"
        if fragment:
            if self.fragment_overlaps is None:
                raise ValueError("Fragment lengths were not set, so there is no fragment level")
            overlaps, max_len = self.fragment_overlaps, self.max_fragment_len
        else:
            overlaps, max_len = self.segment_overlaps, self.max_segment_len

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
        if fragment:
            if self.fragment_overlaps is None:
                raise ValueError("Fragment lengths were not set, so there is no fragment level")
            self.fragment_overlaps._require_candidates().add(starts, ends, extra_penalty)
        else:
            self.segment_overlaps._require_candidates().add(starts, ends, extra_penalty)

    def chop(self) -> Assembly:
        """Cut the sequence into the fewest segments whose overlaps all satisfy the constraints."""

        if self.segment_overlaps.candidates is None:
            self.generate_candidates(fragment=False)
        self.gap_check(fragment=False)
        if self.fragment_overlaps is not None:
            if self.fragment_overlaps.candidates is None:
                self.generate_candidates(fragment=True)
            self.gap_check(fragment=True)

        c = self.segment_overlaps.candidates
        starts, ends = c.starts, c.ends

        weights, back = self._shortest_path(starts, ends, self.min_segment_len, self.max_segment_len)
        if not np.isfinite(weights[-1]):
            raise NoAssemblyError(self._why_stuck(starts, ends, weights))

        path = self._walk_back(back)
        return Assembly(
            segments=[(int(starts[i]), int(ends[j])) for i, j in zip(path, path[1:])],
            overlaps=[(int(starts[i]), int(ends[i])) for i in path[1:-1]],
        )

    @staticmethod
    def _shortest_path(starts: np.ndarray, ends: np.ndarray,
                       min_len: int, max_len: int) -> tuple[np.ndarray, np.ndarray]:
        """Fewest-edges path over the overlaps, each edge one segment.

        The overlaps are sorted by start and every edge runs forward, so one sweep in that
        order settles each vertex before it is used. An edge from i to j is the segment
        [start_i, end_j); it carries weight 1, so a shortest path is a shortest assembly.

        Each vertex relaxes its whole window at once. Starts are sorted, so the window ends
        at the first start out of segment reach, and the segment length is then checked
        against both bounds per candidate.
        """
        n = len(starts)
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
            weight = weights[i] + 1.0
            better = (length >= min_len) & (length <= max_len) & (weight < weights[window])

            relaxed = np.flatnonzero(better) + i + 1
            weights[relaxed] = weight
            back[relaxed] = i

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
