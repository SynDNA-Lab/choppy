from dataclasses import dataclass, field, InitVar
import numpy as np
from choppy.constraints import AbstractCondition


class NoAssemblyError(RuntimeError):
    """Raised when no chain of overlaps spans the sequence."""


@dataclass
class Overlaps:
    """Every potential overlap of a sequence, together with the conditions scored on it.

    Each condition's map carries the grid that converts between an overlap (start, end)
    and the cell that scores it, so a candidate overlap is a cell. Candidates are the
    cells that satisfy all constraints at once.
    """
    seq: str
    min_len: int # miminal overlap length
    max_len: int # maximal overlap length
    min_step: int | None = None # thinning parameter for overlaps (delta_start + delta_end >= min_step)

    objectives: list[AbstractCondition] = field(default_factory=list)
    constraints: list[AbstractCondition] = field(default_factory=list)
    candidates: list[tuple[int, int]] = field(default_factory=list)

    shape: tuple[int, int] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self.shape = (len(self.seq) - self.min_len + 1, self.max_len - self.min_len + 1)
        if self.min_step is not None and self.min_step <= 0:
            raise ValueError(f"min_step must be positive, got {self.min_step}")
        if self.min_step == 1:
            self.min_step = None  # no thinning needed if min_step is 1

    def add_condition(self, condition: AbstractCondition) -> None:
        condition.build_map(self.seq, self.min_len, self.max_len)
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
        for i in range(self.min_step // 2):
            for j in range(self.min_step - i * 2):
                offsets.append((i, j))
        offsets.pop(0)  # don't include (0, 0) because it would always be kept
        offsets = np.array(offsets)

        for i in range(mask.shape[0]):
            for j in range(mask.shape[1]):
                if not mask[i, j]:
                    continue

                rows, cols = i + offsets[:, 0], j + offsets[:, 1]
                inside = (rows < mask.shape[0]) & (cols < mask.shape[1])
                mask[rows[inside], cols[inside]] = False
        return mask

    def generate_candidates(self) -> list[tuple[int, int]]:
        """Potential overlaps as (start, end), sorted by start then end."""
        mask = self._thin(self.feasible_mask())

        rows, cols = np.nonzero(mask)
        self.candidates = list(zip(rows.tolist(), (rows + self.min_len + cols).tolist()))

        return self.candidates


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
    min_segment_len: InitVar[int]
    max_segment_len: InitVar[int]
    min_fragment_overlap_len: InitVar[int | None] = None
    max_fragment_overlap_len: InitVar[int | None] = None
    min_fragment_len: InitVar[int | None] = None
    max_fragment_len: InitVar[int | None] = None
    min_step: int | None = None

    fragments: Overlaps | None = field(init=False, repr=False)
    segments: Overlaps = field(init=False, repr=False)

    def __post_init__(self, min_segment_overlap_len, max_segment_overlap_len, min_segment_len, max_segment_len,
                 min_fragment_overlap_len, max_fragment_overlap_len, min_fragment_len, max_fragment_len):
        self.seq = self.seq.upper()
        self.seq_len = len(self.seq)

        self._check_lengths(
            "segment",
            min_segment_overlap_len,
            max_segment_overlap_len,
            min_segment_len,
            max_segment_len,
        )

        fragment_params = (
            min_fragment_overlap_len,
            max_fragment_overlap_len,
            min_fragment_len,
            max_fragment_len,
        )
        n_set = sum(p is not None for p in fragment_params)
        if n_set not in (0, 4):
            raise ValueError("Fragment lengths must be either all set or all None")
        self.use_fragments = n_set == 4
        if self.use_fragments:
            self._check_lengths("fragment", *fragment_params)

        self.min_segment_len = min_segment_len
        self.max_segment_len = max_segment_len
        self.min_fragment_len = min_fragment_len
        self.max_fragment_len = max_fragment_len

        self.fragments = Overlaps(self.seq, min_fragment_overlap_len, max_fragment_overlap_len,
                                  self.min_step) if self.use_fragments else None
        self.segments = Overlaps(self.seq, min_segment_overlap_len, max_segment_overlap_len,
                                 self.min_step)

    @staticmethod
    def _check_lengths(name: str, min_overlap_len, max_overlap_len, min_len, max_len) -> None:
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
        if min_len <= 2 * min_overlap_len:
            raise ValueError(
                f"min_{name}_len ({min_len}) must be > 2 * min_{name}_overlap_len ({min_overlap_len})"
            )
        if max_len <= 2 * max_overlap_len:
            raise ValueError(
                f"max_{name}_len ({max_len}) must be > 2 * max_{name}_overlap_len ({max_overlap_len})"
            )
        if max_len <= min_len:
            raise ValueError(f"max_{name}_len ({max_len}) must be > min_{name}_len ({min_len})")

    def add_condition(self, condition: AbstractCondition, fragment: bool = False) -> None:
        if fragment:
            if self.fragments is None:
                raise ValueError("Fragment lengths were not set, so there is no fragment level")
            self.fragments.add_condition(condition)
        else:
            self.segments.add_condition(condition)

    def cut(self) -> Assembly:
        """Cut the sequence into the fewest segments whose overlaps all satisfy the constraints."""
        self.segments.generate_candidates()
        starts, ends = self._with_sentinels(self.segments)

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
