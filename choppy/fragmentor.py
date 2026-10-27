from dataclasses import dataclass, field, InitVar
from choppy.constraints import AbstractCondition

@dataclass
class Overlaps:
    seq: str
    min_len: int
    max_len: int

    objectives: list[AbstractCondition] = field(default_factory=list)
    constraints: list[AbstractCondition] = field(default_factory=list)
    candidates: list[tuple[int, int]] = field(default_factory=list)

    def add_condition(self, condition: AbstractCondition) -> None:
        condition.build_map(self.seq, self.min_overlap_len, self.max_overlap_len)
        if condition.opt_val is not None:
            self.objectives.append(condition)
        if condition.min_val is not None or condition.max_val is not None:
            self.constraints.append(condition)

    def generate_candidates(self, seq:str) -> None:
        self.candidates = [(start, length) for start in range(len(seq) - self.min_len + 1)
                           for length in range(self.min_len, self.max_len + 1)
                           if start + length <= len(seq)]

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

    fragments: Overlaps = field(init=False, repr=False)
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

        self.fragments = Overlaps(self.seq, min_fragment_overlap_len, max_fragment_overlap_len)
        self.segments = Overlaps(self.seq, min_segment_overlap_len, max_segment_overlap_len)

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

    def add_condition(self, condition: AbstractCondition, fragment: bool) -> None:
        if fragment:
            self.fragments.add_condition(condition)
        else:
            self.segments.add_condition(condition)

