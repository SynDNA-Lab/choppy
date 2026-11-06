from dataclasses import dataclass, field
from abc import ABC, abstractmethod
import numpy as np
from Bio import SeqIO
import re

def reverse_complement(seq: str) -> str:
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]

def encode_kmers(seq: str, kmer: int) -> np.ndarray:
    lut = np.full(256, 255, np.uint8)
    for i, b in enumerate(b"ACGT"): lut[b] = i
    codes = lut[np.frombuffer(seq.encode(), np.uint8)]
    m = len(codes) - kmer + 1
    enc = np.zeros(m, dtype=np.uint64)
    for j in range(kmer):
        enc = (enc << np.uint64(2)) | codes[j:j + m]
    return enc

@dataclass(frozen=True)
class OverlapGrid:
    """Converts between overlap coordinates and condition map cells.

    Row i is the overlap start, column j the overlap length offset from min_len, so cell
    (i, j) is the overlap (i, i + min_len + j). min_len and max_len are overlap lengths,
    not segment or fragment lengths.

    `cell` and `overlap` take scalars or arrays. The whole grid at once is `starts` (a
    column), `lengths` (a row) and `ends` (the full grid the two broadcast to), so a rule
    over every overlap is written as one expression instead of rebuilt index arithmetic.
    """
    seq_len: int
    min_len: int
    max_len: int

    @property
    def shape(self) -> tuple[int, int]:
        return (self.seq_len - self.min_len + 1, self.max_len - self.min_len + 1)

    def cell(self, start, end):
        """Convert an overlap (start, end) to a cell in the condition map."""
        length = np.subtract(end, start)
        if np.any(length < self.min_len) or np.any(length > self.max_len):
            raise ValueError(f"Overlap length {length} is out of bounds [{self.min_len}, {self.max_len}]")
        i = start
        j = length - self.min_len
        return i, j

    def overlap(self, i, j):
        """Convert a cell in the condition map to an overlap (start, end)."""
        n_starts, n_lens = self.shape
        if np.any(np.less(i, 0)) or np.any(np.greater_equal(i, n_starts)) \
                or np.any(np.less(j, 0)) or np.any(np.greater_equal(j, n_lens)):
            raise ValueError(f"Cell ({i}, {j}) is out of bounds {self.shape}")
        start = i
        end = start + self.min_len + j
        return start, end

    @property
    def starts(self) -> np.ndarray:
        """Overlap start of each row, shaped (n_starts, 1) to broadcast over the grid."""
        return np.arange(self.shape[0], dtype=np.int64)[:, None]

    @property
    def lengths(self) -> np.ndarray:
        """Overlap length of each column, shaped (1, n_lens) to broadcast over the grid."""
        return np.arange(self.min_len, self.max_len + 1, dtype=np.int64)[None, :]

    @property
    def ends(self) -> np.ndarray:
        """Overlap end of every cell, shaped (n_starts, n_lens)."""
        return self.starts + self.lengths

    @property
    def inside(self) -> np.ndarray:
        """Cells whose overlap ends within the sequence."""
        return self.ends <= self.seq_len

    def edge_mask(self, min_edge_len: int) -> np.ndarray:
        """Cells leaving at least min_edge_len of sequence beyond the overlap at both ends."""
        return (self.ends >= min_edge_len) & (self.starts <= self.seq_len - min_edge_len)

@dataclass
class OverlapMap:
    grid: OverlapGrid
    values: np.ndarray

    def get_cell(self, i: int, j: int) -> np.float16:
        return self.values[i, j]

    def get_overlap(self, start: int, end: int) -> np.float16:
        i, j = self.grid.cell(start, end)
        return self.values[i, j]

@dataclass
class AbstractCondition(ABC):
    min_val: np.float16 | None = None
    max_val: np.float16 | None = None
    opt_val: np.float16 | None = None
    grid: OverlapGrid | None = field(default=None, init=False, repr=False)
    map: OverlapMap | None = field(default=None, init=False, repr=False)

    def bind(self, grid: OverlapGrid) -> None:
        """Attach the grid this condition is scored on, normally by the owning Overlaps."""
        if self.grid is not None and self.grid != grid:
            raise ValueError(
                f"Condition is already bound to {self.grid}, cannot rebind to {grid}. "
                f"Use a separate instance for the other level."
            )
        self.grid = grid

    def build_map(self, seq: str, grid: OverlapGrid | None = None) -> OverlapMap:
        """Score this condition over every overlap of the grid it is bound to."""
        if grid is not None:
            self.bind(grid)
        if self.grid is None:
            raise ValueError("No grid bound: pass one to build_map, or call bind first")
        if len(seq) != self.grid.seq_len:
            raise ValueError(
                f"Sequence is {len(seq)} long but the bound grid is for {self.grid.seq_len}"
            )
        self.map = OverlapMap(self.grid, np.full(self.grid.shape, 0.0, dtype=np.float16))
        self._fill_map(seq)
        return self.map

    @abstractmethod
    def _fill_map(self, seq: str) -> None: ...

    def feasible(self) -> np.ndarray:
        """Boolean mask of the overlaps that satisfy this condition as a constraint.

        min_val is a floor, max_val a ceiling; a condition with neither set is not a
        constraint and admits everything. NaN entries (overlaps running past the end of
        the sequence) are never feasible.
        """
        if self.map is None:
            raise ValueError("build_map must be called before feasible")
        values = self.map.values
        ok = ~np.isnan(values)
        if self.min_val is not None:
            ok &= values >= self.min_val
        if self.max_val is not None:
            ok &= values <= self.max_val
        return ok

## TODO: rewrite against the condition's own thresholds (max_val) and fill the whole map
# class PolyX(AbstractCondition):
#
#     def score(self, seq: str) -> np.float16:
#         return np.float16(len(seq) > self.max_len)
#
#     def _fill_map(self, seq: str) -> None:
#         self.map.values[0, 0] = self.score(seq)

@dataclass
class BackgroundKmers:
    kmer: int
    keys: np.ndarray
    counts: np.ndarray

    def __post_init__(self):
        self._check_kmer(self.kmer)

    @staticmethod
    def _check_kmer(kmer: int) -> None:
        if kmer < 2 or kmer > 16:
            raise ValueError(f"kmer must be between 2 and 16, got {kmer}")


    @classmethod
    def from_fasta(cls, fasta: str, kmer: int) -> "BackgroundKmers":
        cls._check_kmer(kmer)
        sequences = [str(record.seq).upper() for record in SeqIO.parse(fasta, "fasta")]
        keys, counts = cls.encode_seqs(sequences, kmer)
        return cls(kmer=kmer, keys=keys, counts=counts)

    @classmethod
    def load(cls, path: str) -> "BackgroundKmers":
        data = np.load(path, allow_pickle=True)
        cls._check_kmer(data["kmer"])
        return cls(kmer=data["kmer"], keys=data["keys"], counts=data["counts"])

    @classmethod
    def from_seqs(cls, sequences: list[str], kmer: int) -> "BackgroundKmers":
        cls._check_kmer(kmer)
        keys, counts = cls.encode_seqs(sequences, kmer)
        return cls(kmer=kmer, keys=keys, counts=counts)

    def save(self, path: str) -> None:
        np.savez_compressed(path, kmer=self.kmer, keys=self.keys, counts=self.counts)

    @staticmethod
    def encode_seqs(sequences: list[str] | str, kmer: int) -> np.ndarray:
        if isinstance(sequences, str):
            sequences = [sequences]

        sequences = sequences + [reverse_complement(s) for s in sequences]

        all_codes = []
        for seq in sequences:
            enc = encode_kmers(seq, kmer)
            all_codes.append(enc)

        keys, counts = np.unique(np.concatenate(all_codes), return_counts=True)
        return keys, counts

    def lookup(self, seq: str) -> np.ndarray:
        enc_seq = encode_kmers(seq, self.kmer)
        idx = np.searchsorted(self.keys, enc_seq)
        idx = np.minimum(idx, len(self.keys) - 1)
        found = self.keys[idx] == enc_seq
        counts = np.zeros(len(enc_seq), dtype=np.int32)
        counts[found] = self.counts[idx[found]]
        return counts

@dataclass
class KmerHomology(AbstractCondition):
    kmer: int = 15

    @abstractmethod
    def _kmer_scores(self, seq: str) -> np.ndarray: ...

    def _fill_map(self, seq: str) -> None:
        grid = self.grid
        scores = self._kmer_scores(seq)
        padded = np.concatenate([scores, np.full(grid.max_len - grid.min_len, np.nan, dtype=np.float16)])
        windows = np.lib.stride_tricks.sliding_window_view(padded, grid.max_len - self.kmer + 1)
        self.map.values[:] = np.maximum.accumulate(windows, axis=1)[:, grid.min_len - self.kmer:]

@dataclass
class BGHomology(KmerHomology):
    bg_fasta: str | None = None
    bg_npz: str | None = None
    bg: BackgroundKmers | None = field(default=None, init=False, repr=False)
    multiplier: float = field(init=False, repr=False)

    def __post_init__(self):
        self._set_bg()
        self.multiplier = np.max(self.bg.counts)
        if self.kmer != self.bg.kmer:
            raise ValueError(f"kmer mismatch: {self.kmer} != {self.bg.kmer}")

    def _set_bg(self) -> BackgroundKmers:
        if self.bg_fasta is None:
            if self.bg_npz is None:
                if self.bg is None:
                    raise ValueError("No background provided")
            else:
                self.bg = BackgroundKmers.load(self.bg_npz)
        else:
            self.bg = BackgroundKmers.from_fasta(self.bg_fasta, self.kmer)

    def _kmer_scores(self, seq: str) -> np.ndarray:
        return np.float16(self.bg.lookup(seq) / self.multiplier)

@dataclass
class SelfHomology(KmerHomology):

    def _kmer_scores(self, seq: str) -> np.ndarray:
        repeats = BackgroundKmers.from_seqs([seq], self.kmer).lookup(seq) - 1
        return np.float16(repeats / max(repeats.max(), 1))

@dataclass
class BoundaryMotif(AbstractCondition):
    motif: str = "GC"

    def _fill_map(self, seq: str) -> None:
        grid = self.grid
        n_starts = grid.shape[0]
        matches = list(re.finditer(self.motif, seq))
        starts = np.array([m.start() for m in matches], dtype=np.int64)
        ends = np.array([m.end() for m in matches], dtype=np.int64)

        opens = np.zeros(n_starts, dtype=bool)
        opens[starts[starts < n_starts]] = True

        closes = np.zeros(grid.seq_len + 1, dtype=bool)
        closes[ends] = True

        stop = grid.ends
        inside = grid.inside
        shuts = closes[np.minimum(stop, grid.seq_len)]

        bounded = opens[:, None] & shuts
        self.map.values[:] = np.where(inside, ~bounded, np.nan).astype(np.float16)

