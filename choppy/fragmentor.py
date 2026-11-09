from dataclasses import dataclass, field, InitVar
import warnings
import numpy as np
from choppy.constraints import AbstractCondition, OverlapGrid, check_seq


class NoAssemblyError(RuntimeError):
    """Raised when no chain of overlaps spans the sequence."""


@dataclass
class Overlaps:
    """Every potential overlap of a sequence, together with the conditions scored on it."""
    seq: str
    min_len: int # miminal overlap length
    max_len: int # maximal overlap length
    min_step: int | None = None # thinning parameter for overlaps (delta_start + delta_end >= min_step)

    objectives: list[AbstractCondition] = field(default_factory=list)
    constraints: list[AbstractCondition] = field(default_factory=list)
    candidates: list[tuple[int, int]] = field(default_factory=list)

    grid: OverlapGrid = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self.grid = OverlapGrid(len(self.seq), self.min_len, self.max_len)
        if self.min_step == 1:
            self.min_step = None  # no thinning needed if min_step is 1

    @property
    def shape(self) -> tuple[int, int]:
        return self.grid.shape

    def add_condition(self, condition: AbstractCondition) -> None:
        condition.build_map(self.seq, self.grid)
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

    def generate_candidates(self, min_edge_len: int) -> list[tuple[int, int]]:
        """Potential overlaps as (start, end), sorted by start then end."""
        mask = self.feasible_mask() & self.grid.edge_mask(min_edge_len)
        mask = self._thin(mask)

        starts, ends = self.grid.overlap(*np.nonzero(mask))
        self.candidates = list(zip(starts.tolist(), ends.tolist()))

        return self.candidates

    def gap_lengths(self) -> np.ndarray:
        """The lengths of the potential consecutive edges."""
        if not self.candidates:
            raise ValueError("No candidates have been generated yet")
        starts, ends = np.array(self.candidates, dtype=np.int64).T
        return ends[1:] - starts[:-1]

    def plot_gap_lengths(self, min_gap_len: int = None, ax=None):
        """Histogram of the lengths of the gaps between consecutive overlaps."""
        if not self.candidates:
            raise ValueError("No candidates have been generated yet")
        lengths = self.gap_lengths()
        if min_gap_len is not None:
            lengths = lengths[lengths >= min_gap_len]
        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
        ax.hist(lengths, bins=range(1, lengths.max() + 2), align="left", color="tab:blue")
        ax.set_xlabel("gap length")
        ax.set_ylabel("count")
        ax.set_title(f"gaps between {len(self.candidates)} candidate overlaps")
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
            ax.set_title(f"[{i}] {type(constraint).__name__}")
            ax.set_xlabel("")
        axes[-1].set_xlabel("overlap start")
        fig.suptitle(f"Overlaps starting in [{a}, {b})")
        fig.tight_layout()
        return axes

    def get_gapped_regions(self, min_gap_len: int) -> list[tuple[int, int]]:
        """The regions of the sequence that are not covered by any candidate overlap.

        Returns a list of (start, end) tuples, sorted by start. Each region is at least
        min_gap_len long.
        """
        if not self.candidates:
            raise ValueError("No candidates have been generated yet")
        starts, ends = np.array(self.candidates, dtype=np.int64).T
        gaps = []
        prev_end = 0
        for start, end in zip(starts, ends):
            if start - prev_end >= min_gap_len:
                gaps.append((prev_end, start))
            prev_end = max(prev_end, end)
        if len(self.seq) - prev_end >= min_gap_len:
            gaps.append((prev_end, len(self.seq)))
        return gaps
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

        self.fragment_overlaps = Overlaps(self.seq, min_fragment_overlap_len, max_fragment_overlap_len,
                                          self.min_step) if self.use_fragments else None
        self.segment_overlaps = Overlaps(self.seq, min_segment_overlap_len, max_segment_overlap_len,
                                         self.min_step)

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

    def add_condition(self, condition: AbstractCondition, fragment: bool = False) -> None:
        if fragment:
            if self.fragment_overlaps is None:
                raise ValueError("Fragment lengths were not set, so there is no fragment level")
            self.fragment_overlaps.add_condition(condition)
        else:
            self.segment_overlaps.add_condition(condition)

    def generate_candidates(self, fragment: bool = False) -> list[tuple[int, int]]:
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

    def chop(self) -> Assembly:
        """Cut the sequence into the fewest segments whose overlaps all satisfy the constraints."""
        
        self.generate_candidates(fragment=False)
        if self.fragment_overlaps is not None:
            self.generate_candidates(fragment=True)

        starts, ends = self._with_sentinels(self.segment_overlaps)

        weights, back = self._shortest_path(starts, ends, self.min_segment_len, self.max_segment_len)
        if not np.isfinite(weights[-1]):
            raise NoAssemblyError(self._why_stuck(starts, ends, weights))

        path = self._walk_back(back)
        return Assembly(
            segments=[(int(starts[i]), int(ends[j])) for i, j in zip(path, path[1:])],
            overlaps=[(int(starts[i]), int(ends[i])) for i in path[1:-1]],
        )

    def _with_sentinels(self, overlaps: Overlaps) -> tuple[np.ndarray, np.ndarray]:
        """Candidates bracketed by zero-length overlaps at each end of the sequence."""
        candidates = np.array(overlaps.candidates, dtype=np.int64).reshape(-1, 2)
        starts = np.concatenate(([0], candidates[:, 0], [self.seq_len]))
        ends = np.concatenate(([0], candidates[:, 1], [self.seq_len]))
        return starts, ends

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
