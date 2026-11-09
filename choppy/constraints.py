from dataclasses import dataclass, field
from abc import ABC, abstractmethod
import warnings
import numpy as np
from Bio import SeqIO
import re

ALPHABET = "ACGT"

def reverse_complement(seq: str) -> str:
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]

def invalid_bases(seq: str) -> np.ndarray:
    """Positions of the characters of seq that are not A, C, G or T."""
    known = np.zeros(256, dtype=bool)
    for b in ALPHABET.encode(): known[b] = True
    return np.flatnonzero(~known[np.frombuffer(seq.encode(), np.uint8)])

def check_seq(seq: str, name: str = "Sequence") -> str:
    """Return seq unchanged, or raise if it holds anything but A, C, G and T.

    Ambiguity codes are not supported anywhere in choppy: `encode_kmers` has only two bits
    per base, so an unknown base would silently alias a real kmer. Callers that take a
    sequence from outside uppercase it first.
    """
    bad = invalid_bases(seq)
    if bad.size:
        found = sorted({seq[p] for p in bad})
        raise ValueError(
            f"{name} contains {bad.size} character(s) outside {ALPHABET}: {found}, "
            f"first at position {int(bad[0])}. Ambiguity codes are not supported."
        )
    return seq

def valid_kmers(seq: str, kmer: int) -> np.ndarray:
    """Mask of the kmers of seq, True where the kmer is made only of A, C, G and T."""
    bad = np.zeros(len(seq), dtype=np.int64)
    bad[invalid_bases(seq)] = 1
    cumulative = np.concatenate(([0], np.cumsum(bad)))
    return (cumulative[kmer:] - cumulative[:-kmer]) == 0

def encode_kmers(seq: str, kmer: int) -> np.ndarray:
    """Base-4 code of every kmer of seq.

    Only A, C, G and T are encoded; any other character makes the codes of the kmers
    covering it meaningless, so callers either reject such a sequence with `check_seq` or
    drop those kmers with `valid_kmers`.
    """
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

    def check_cell(self, i, j) -> None:
        """Raise unless (i, j) is a cell of this grid. Takes scalars or arrays."""
        n_starts, n_lens = self.shape
        if np.any(np.less(i, 0)) or np.any(np.greater_equal(i, n_starts)) \
                or np.any(np.less(j, 0)) or np.any(np.greater_equal(j, n_lens)):
            raise ValueError(f"Cell ({i}, {j}) is out of bounds {self.shape}")

    def cell(self, start, end):
        """Convert an overlap (start, end) to a cell in the condition map."""
        length = np.subtract(end, start)
        if np.any(length < self.min_len) or np.any(length > self.max_len):
            raise ValueError(f"Overlap length {length} is out of bounds [{self.min_len}, {self.max_len}]")
        i = start
        j = length - self.min_len
        # the end is deliberately not bounded: the grid holds cells whose overlap runs past
        # the end of the sequence, and edge_mask is what drops them.
        self.check_cell(i, j)
        return i, j

    def overlap(self, i, j):
        """Convert a cell in the condition map to an overlap (start, end)."""
        self.check_cell(i, j)
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

    def edge_mask(self, min_edge_len: int) -> np.ndarray:
        """
        Cells that represent overlaps that can never potentially be used since they
        correspond to a segment/fragment that is shroter than the provided minimum.
        """
        return (self.ends >= min_edge_len) & (self.starts <= self.seq_len - min_edge_len)

@dataclass
class OverlapMap:
    grid: OverlapGrid
    values: np.ndarray

    def get_cell(self, i: int, j: int) -> np.float16:
        self.grid.check_cell(i, j)
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
        """Score this condition over every overlap of the grid it is bound to.

        Scores run from 0 to 1 (times an optional multiplier). 0 is always the best score,
        1 is the worst.
        """
        if self.min_val is None and self.max_val is None and self.opt_val is None:
            raise ValueError("Condition has no thresholds or objectives set: it is neither a constraint nor an objective "
                             "and build_map is meaningless")
        if grid is not None:
            self.bind(grid)
        if self.grid is None:
            raise ValueError("No grid bound: pass one to build_map, or call bind first")
        if len(seq) != self.grid.seq_len:
            raise ValueError(
                f"Sequence is {len(seq)} long but the bound grid is for {self.grid.seq_len}"
            )
        self.map = OverlapMap(self.grid, np.full(self.grid.shape, 1.0, dtype=np.float16))
        self._fill_map(seq)
        return self.map

    @abstractmethod
    def _fill_map(self, seq: str) -> None: ...

    def feasible(self) -> np.ndarray:
        """Boolean mask of the overlaps that satisfy this condition as a constraint.

        min_val is a floor, max_val a ceiling; a condition with neither set is not a
        constraint and admits everything. A NaN score fails both comparisons, so it is
        infeasible without a separate check for it.
        """
        if self.map is None:
            raise ValueError("build_map must be called before feasible")
        values = self.map.values
        if self.min_val is None and self.max_val is None:
            return np.ones(values.shape, dtype=bool)
        if self.min_val is None:
            return values <= self.max_val
        ok = values >= self.min_val
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
        """Count every kmer of the sequences and of their reverse complements.

        Unlike the sequence being cut, a background may legitimately hold ambiguity codes
        and assembly gaps. They cannot be encoded, so the kmers covering them are dropped
        and reported; the rest of the background is still counted.
        """
        if isinstance(sequences, str):
            sequences = [sequences]
        sequences = [s.upper() for s in sequences]

        too_short = [len(s) for s in sequences if len(s) < kmer]
        if too_short:
            warnings.warn(
                f"{len(too_short)} background sequence(s) are shorter than kmer {kmer} "
                f"and hold no kmer at all; they are ignored",
                stacklevel=2,
            )
        sequences = [s for s in sequences if len(s) >= kmer]
        if not sequences:
            raise ValueError(f"No background sequence is as long as kmer {kmer}")

        all_codes = []
        n_dropped = 0
        for seq in sequences:
            for strand in (seq, reverse_complement(seq)):
                ok = valid_kmers(strand, kmer)
                all_codes.append(encode_kmers(strand, kmer)[ok])
            n_dropped += int((~valid_kmers(seq, kmer)).sum())

        if n_dropped:
            n_bases = sum(invalid_bases(s).size for s in sequences)
            warnings.warn(
                f"Background holds {n_bases} base(s) outside {ALPHABET}; the "
                f"{n_dropped} kmer(s) covering them are ignored (per strand)",
                stacklevel=2,
            )

        pooled = np.concatenate(all_codes)
        if pooled.size == 0:
            raise ValueError(f"No valid {kmer}-mer in the background: every kmer holds a base outside {ALPHABET}")
        keys, counts = np.unique(pooled, return_counts=True)
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
    cap: int = 2048

    @abstractmethod
    def _kmer_scores(self, seq: str) -> np.ndarray: ...

    def _fill_map(self, seq: str) -> None:
        """Worst kmer score over each overlap, as a running max along the length axis.

        The pad gives the sliding window one row per overlap start; it is filled with 1.0,
        the worst score, and the running max reaches it exactly on the cells whose overlap
        ends past the end of the sequence.

        An overlap shorter than kmer holds no complete kmer, so there is no evidence of
        homology against it and it keeps the best score, 0. Those are the first
        `kmer - min_len` length columns, which have no window to slide over.
        """
        n_cols = self.grid.shape[1]
        n_short = min(max(self.kmer - self.grid.min_len, 0), n_cols)
        self.map.values[:, :n_short] = 0.0
        if n_short == n_cols:  # kmer > max_len: no overlap is long enough to hold one
            return

        scores = self._kmer_scores(seq)
        padded = np.concatenate([scores, np.full(self.grid.max_len - self.grid.min_len, 1.0, dtype=np.float16)])
        windows = np.lib.stride_tricks.sliding_window_view(padded, self.grid.max_len - self.kmer + 1)
        accumulated = np.maximum.accumulate(windows, axis=1)
        # column j is the overlap of length min_len + j, whose last kmer starts at
        # offset min_len + j - kmer; the short columns above are the negative ones.
        self.map.values[:, n_short:] = accumulated[:, max(self.grid.min_len - self.kmer, 0):]

@dataclass
class BGHomology(KmerHomology):
    bg_fasta: str | None = None
    bg_npz: str | None = None
    bg: BackgroundKmers | None = field(default=None, init=False, repr=False)

    def __post_init__(self):
        if self.cap < 1 or self.cap > 2048:
            raise ValueError(f"cap must be between 1 and 2048, got {self.cap}")
        self._set_bg()
        if self.kmer != self.bg.kmer:
            raise ValueError(f"kmer mismatch: {self.kmer} != {self.bg.kmer}")
        max_counts = np.max(self.bg.counts)
        if max_counts < self.cap:
            self.cap = int(max_counts)
        else:
            self.bg.counts = np.minimum(self.bg.counts, self.cap)

    def _set_bg(self) -> None:
        """Load the background: from the FASTA if there is one, otherwise from the npz."""
        if self.bg_fasta is None and self.bg_npz is None:
            raise ValueError("No background provided: set either bg_fasta or bg_npz")
        if self.bg_fasta is not None:
            if self.bg_npz is not None:
                warnings.warn(
                    "both bg_fasta and bg_npz are set; bg_fasta takes priority and "
                    f"{self.bg_npz} is ignored",
                    stacklevel=3,
                )
            self.bg = BackgroundKmers.from_fasta(self.bg_fasta, self.kmer)
        else:
            self.bg = BackgroundKmers.load(self.bg_npz)

    def _kmer_scores(self, seq: str) -> np.ndarray:
        return np.float16(self.bg.lookup(seq) / self.cap)

@dataclass
class SelfHomology(KmerHomology):

    def _kmer_scores(self, seq: str) -> np.ndarray:
        repeats = BackgroundKmers.from_seqs([seq], self.kmer).lookup(seq) - 1
        return np.float16(repeats / max(repeats.max(), 1))

@dataclass
class BoundaryMotif(AbstractCondition):
    """Overlaps that both open and close on a match of the motif score 0, the rest 1.

    The two ends are independent: the motif opening the overlap and the one closing it are
    separate occurrences. The motif is matched inside a lookahead so that overlapping
    occurrences are all found ("GG" starts at both 0 and 1 of "GGG"); note that only one
    match is taken per position, so a motif that can match several lengths at the same
    place ("[GC]{2,3}") contributes only its greedy end.
    """
    motif: str = "GC"

    def _fill_map(self, seq: str) -> None:
        grid = self.grid
        n_starts = grid.shape[0]
        # the lookahead keeps the match zero-width, so the scan does not skip past an
        # occurrence that starts inside the previous one; group 1 carries the real span.
        matches = list(re.finditer(f"(?=({self.motif}))", seq))
        starts = np.array([m.start(1) for m in matches], dtype=np.int64)
        ends = np.array([m.end(1) for m in matches], dtype=np.int64)

        opens = np.zeros(n_starts, dtype=bool)
        opens[starts[starts < n_starts]] = True

        closes = np.zeros(grid.seq_len + 1, dtype=bool)
        closes[ends] = True

        # clamped because the grid holds cells whose overlap ends past the sequence;
        # their score is meaningless but edge_mask drops them before they are used.
        shuts = closes[np.minimum(grid.ends, grid.seq_len)]

        bounded = opens[:, None] & shuts
        self.map.values[:] = (~bounded).astype(np.float16)

