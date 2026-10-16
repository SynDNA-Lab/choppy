from dataclasses import dataclass, field
from choppy.constraints import AbstractCondition

@dataclass
class Fragmentor:
    seq: str
    min_overlap_len: int
    max_overlap_len: int
    constraints: list[AbstractCondition] = field(default_factory=list)
    objectives: list[AbstractCondition] = field(default_factory=list)

    def __post_init__(self):
        self.seq = self.seq.upper()
        self.seq_len = len(self.seq)

    def add_condition(self, condition: AbstractCondition) -> None:
        condition.build_map(self.seq, self.min_overlap_len, self.max_overlap_len)
        if condition.opt_val is not None:
            self.objectives.append(condition)
        if condition.min_val is not None or condition.max_val is not None:
            self.constraints.append(condition)
