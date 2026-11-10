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

    def slices(self, sub: "OverlapGrid") -> tuple[slice, slice]:
        """The rows and columns of this grid that hold sub, a grid over the same sequence
        whose length range lies inside this one's.

        Cell (i, j) of sub is the overlap (i, i + sub.min_len + j), which is cell
        (i, j + sub.min_len - min_len) here; sub has fewer starts when its min_len is
        larger, and those are the first rows.
        """
        if sub.seq_len != self.seq_len:
            raise ValueError(f"{sub} is over a sequence of {sub.seq_len}, not {self.seq_len}")
        if sub.min_len < self.min_len or sub.max_len > self.max_len:
            raise ValueError(
                f"Overlap lengths [{sub.min_len}, {sub.max_len}] of {sub} are not within "
                f"[{self.min_len}, {self.max_len}]"
            )
        offset = sub.min_len - self.min_len
        return slice(0, sub.shape[0]), slice(offset, offset + sub.shape[1])

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
    multiplier: float = 1.0

    def get_cell(self, i: int, j: int) -> np.float16:
        self.grid.check_cell(i, j)
        return self.values[i, j] * self.multiplier

    def get_overlap(self, start: int, end: int) -> np.float16:
        i, j = self.grid.cell(start, end)
        return self.values[i, j] * self.multiplier

    def get_values(self) -> np.ndarray:
        return self.values * self.multiplier

    def view(self, grid: OverlapGrid) -> "OverlapMap":
        """This map seen through grid, a sub-grid of its own (see OverlapGrid.slices).

        The values are a numpy view, not a copy: the two maps share their memory.
        """
        rows, cols = self.grid.slices(grid)
        return OverlapMap(grid, self.values[rows, cols], self.multiplier)

@dataclass(frozen=True)
class AbstractScorer(ABC):
    """Scores every overlap of a grid: the map-building half of a condition.

    A scorer holds only the parameters that shape its map, never thresholds, so it is
    frozen and hashed by those parameters: two equal scorers would build the same map.
    Subclasses must be frozen dataclasses too.
    """

    def build_map(self, seq: str, grid: OverlapGrid) -> OverlapMap:
        """Score every overlap of grid.

        Scores run from 0 to 1 (times the multiplier _fill returns). 0 is always the best
        score, 1 is the worst.

        Cells whose overlap runs past the end of the sequence always score 1, whatever
        _fill wrote there, so nothing downstream has to filter them out.
        """
        if len(seq) != grid.seq_len:
            raise ValueError(f"Sequence is {len(seq)} long but the grid is for {grid.seq_len}")
        values = np.full(grid.shape, 1.0, dtype=np.float16)
        multiplier = self._fill(seq, grid, values)
        values[grid.ends > grid.seq_len] = 1.0
        return OverlapMap(grid, values, multiplier)

    @abstractmethod
    def _fill(self, seq: str, grid: OverlapGrid, values: np.ndarray) -> float:
        """Write the scores into values (shaped like grid, pre-filled with 1) and return
        the multiplier that turns them back into the scorer's own units."""


THRESHOLDS = ("min_val", "max_val", "opt_val")
_UNSET = object()  # "not given" for Condition.update, where None means "remove"


def _same(a, b) -> bool:
    """Whether two threshold values are the same, None included."""
    if a is None or b is None:
        return a is b
    return a == b


@dataclass(eq=False)
class Condition:
    """A scorer with thresholds: a constraint (min_val and/or max_val), an objective
    (opt_val), or both.

    It compares by identity, unlike its scorer: two conditions on the same scorer are
    still separate conditions. The map is attached by the owning Overlaps. name is
    optional; when set it must be unique on the level, which can then look it up by it.

    The thresholds can be changed after the condition is added, with update or by
    plain assignment (`c.max_val = 3`). Every change that alters a value bumps version,
    which is how the owning level notices that its candidates are out of date. A change
    that would leave no threshold at all is refused and leaves the condition as it was.
    The scorer cannot be swapped once the condition is bound: its map would no longer
    match it.
    """
    scorer: AbstractScorer
    min_val: np.float16 | None = None
    max_val: np.float16 | None = None
    opt_val: np.float16 | None = None
    name: str | None = None
    map: OverlapMap | None = field(default=None, init=False, repr=False)
    version: int = field(default=0, init=False, repr=False)

    def __post_init__(self) -> None:
        self._check_thresholds(self.min_val, self.max_val, self.opt_val)
        # __init__ leaves an init=False field with a plain default to the class attribute,
        # so version only reaches the instance here: its presence in __dict__ marks the
        # end of construction, after which threshold assignments go through update
        object.__setattr__(self, "version", 0)

    @staticmethod
    def _check_thresholds(min_val, max_val, opt_val) -> None:
        if min_val is None and max_val is None and opt_val is None:
            raise ValueError("Condition has no thresholds or objectives set: it is neither a constraint "
                             "nor an objective")

    def __setattr__(self, name: str, value) -> None:
        if "version" not in self.__dict__:  # still inside __init__
            super().__setattr__(name, value)
        elif name in THRESHOLDS:
            self.update(**{name: value})
        elif name == "scorer" and self.map is not None:
            raise AttributeError("The scorer of a bound condition cannot be changed: add a new Condition")
        else:
            super().__setattr__(name, value)

    def update(self, min_val=_UNSET, max_val=_UNSET, opt_val=_UNSET) -> None:
        """Change any of the thresholds at once; None removes one, an omitted one is kept.

        The new set is checked as a whole before anything is changed, so a condition can
        be turned from a constraint into an objective in one call. version is bumped once,
        and only if some value actually changed.
        """
        given = {"min_val": min_val, "max_val": max_val, "opt_val": opt_val}
        new = {k: (getattr(self, k) if v is _UNSET else v) for k, v in given.items()}
        self._check_thresholds(**new)
        changed = {k: v for k, v in new.items() if not _same(v, getattr(self, k))}
        if not changed:
            return
        for k, v in changed.items():
            object.__setattr__(self, k, v)
        object.__setattr__(self, "version", self.version + 1)

    @property
    def grid(self) -> OverlapGrid | None:
        return None if self.map is None else self.map.grid

    @property
    def label(self) -> str:
        """The name if there is one, otherwise the scorer's class name; used in plots."""
        return self.name if self.name is not None else type(self.scorer).__name__

    @property
    def is_constraint(self) -> bool:
        return self.min_val is not None or self.max_val is not None

    @property
    def is_objective(self) -> bool:
        return self.opt_val is not None

    def bind(self, overlap_map: OverlapMap) -> None:
        """Attach the map this condition reads its scores from, normally by the owning
        Overlaps. A condition belongs to one level only, so it is bound once."""
        if self.map is not None:
            raise ValueError(
                f"Condition is already bound to a map on {self.map.grid}. "
                f"Use a separate Condition for the other level."
            )
        self.map = overlap_map

    def feasible(self) -> np.ndarray:
        """Boolean mask of the overlaps that satisfy this condition as a constraint.

        min_val is a floor, max_val a ceiling; a condition with neither set is not a
        constraint and admits everything. A NaN score fails both comparisons, so it is
        infeasible without a separate check for it.
        """
        if self.map is None:
            raise ValueError("Condition must be bound to a map before feasible")
        values = self.map.get_values()
        if self.min_val is None and self.max_val is None:
            return np.ones(values.shape, dtype=bool)
        if self.min_val is None:
            return values <= self.max_val
        ok = values >= self.min_val
        if self.max_val is not None:
            ok &= values <= self.max_val
        return ok

    # ---- visualization -------------------------------------------------------------
    # matplotlib is an optional dependency (the `plot` extra), so it is imported only
    # when a plot is drawn.

    def _plot_axes(self, ax):
        """The axes to draw on: ax itself, or a fresh figure's when ax is None."""
        if self.map is None:
            raise ValueError("Condition must be bound to a map before plotting")
        if ax is None:
            import matplotlib.pyplot as plt
            _, ax = plt.subplots(figsize=(8, 4))
        return ax

    def _plot_thresholds(self, ax) -> None:
        """Horizontal lines at min_val, max_val and opt_val, for those that are set."""
        styles = (
            ("min_val", self.min_val, dict(color="0.45", linestyle="--")),
            ("max_val", self.max_val, dict(color="0.45", linestyle="-.")),
            ("opt_val", self.opt_val, dict(color="tab:green", linestyle=":")),
        )
        for name, value, style in styles:
            if value is not None:
                ax.axhline(float(value), linewidth=1.5, label=f"{name} = {float(value):g}", **style)

    def plot_rows(self, a: int = 0, b: int = None, ax=None):
        """Plot, for each overlap start in [a, b), its minimal score across all overlap lengths.

        Draws on ax if given, otherwise on a new figure, and returns the axes.
        """
        ax = self._plot_axes(ax)
        n_starts = self.grid.shape[0]
        if b is None:
            b = n_starts
        if not 0 <= a < b <= n_starts:
            raise ValueError(f"Row range [{a}, {b}) must be non-empty and within [0, {n_starts})")

        best = self.map.get_values()[a:b].min(axis=1).astype(np.float64)

        ax.scatter(np.arange(a, b), best, s=8, color="tab:blue", label="min score", zorder=3)
        self._plot_thresholds(ax)
        ax.set_xlabel("overlap start")
        ax.set_ylabel("score")
        ax.set_title(f"{self.label}: overlaps starting in [{a}, {b})")
        ax.grid(True, color="0.9")
        ax.legend()
        return ax

class MapStore:
    """One map per scorer over one grid, shared by every condition and level using it.

    A map is built the first time its scorer is asked for; equal scorers (same type and
    parameters) get the same map. The stored values are read-only, since a write through
    one condition would silently change every other condition reading the same map.
    """

    def __init__(self, seq: str, grid: OverlapGrid):
        if len(seq) != grid.seq_len:
            raise ValueError(f"Sequence is {len(seq)} long but the grid is for {grid.seq_len}")
        self.seq = seq
        self.grid = grid
        self._maps: dict[AbstractScorer, OverlapMap] = {}

    def get(self, scorer: AbstractScorer) -> OverlapMap:
        """The map of scorer over the whole grid, built on first use."""
        overlap_map = self._maps.get(scorer)
        if overlap_map is None:
            overlap_map = scorer.build_map(self.seq, self.grid)
            overlap_map.values.flags.writeable = False
            self._maps[scorer] = overlap_map
        return overlap_map

    def view(self, scorer: AbstractScorer, grid: OverlapGrid) -> OverlapMap:
        """The map of scorer seen through grid, a sub-grid of the store's."""
        return self.get(scorer).view(grid)

    def __contains__(self, scorer: AbstractScorer) -> bool:
        return scorer in self._maps

    def __len__(self) -> int:
        return len(self._maps)

    def __repr__(self) -> str:
        return f"MapStore({len(self)} maps on {self.grid})"

## TODO: rewrite as a scorer (thresholds live on the Condition) and fill the whole map
# class PolyX(AbstractScorer):
#
#     def score(self, seq: str) -> np.float16:
#         return np.float16(len(seq) > self.max_len)
#
#     def _fill(self, seq, grid, values) -> float:
#         values[0, 0] = self.score(seq)
#         return 1.0

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

@dataclass(frozen=True)
class KmerHomology(AbstractScorer):
    kmer: int = 15

    @abstractmethod
    def _kmer_scores(self, seq: str) -> tuple[np.ndarray, float]: ...

    def _fill(self, seq: str, grid: OverlapGrid, values: np.ndarray) -> float:
        """Worst kmer score over each overlap, as a running max along the length axis.

        An overlap shorter than kmer holds no complete kmer, so there is no evidence of
        homology against it and it keeps the best score, 0. Those are the first
        `kmer - min_len` length columns, which have no window to slide over.

        The kmer scores are computed even when no overlap is long enough to need them:
        the multiplier comes with them, and it belongs to the scorer and the sequence, not
        to the grid, so a map and any view of a wider map must agree on it.
        """
        n_cols = grid.shape[1]
        n_short = min(max(self.kmer - grid.min_len, 0), n_cols)
        values[:, :n_short] = 0.0
        scores, multiplier = self._kmer_scores(seq)
        if n_short == n_cols:  # kmer > max_len: no overlap is long enough to hold one
            return multiplier

        padded = np.concatenate([scores, np.full(grid.max_len - grid.min_len, 1.0, dtype=np.float16)])
        windows = np.lib.stride_tricks.sliding_window_view(padded, grid.max_len - self.kmer + 1)
        accumulated = np.maximum.accumulate(windows, axis=1)
        # column j is the overlap of length min_len + j, whose last kmer starts at
        # offset min_len + j - kmer; the short columns above are the negative ones.
        values[:, n_short:] = accumulated[:, max(grid.min_len - self.kmer, 0):]
        return multiplier

@dataclass(frozen=True)
class BGHomology(KmerHomology):
    """Kmer homology against a background (e.g. the host genome), counts capped at cap.

    The loaded background is not part of equality or the hash: two scorers on the same
    background source, kmer and cap are equal. cap is lowered to the highest count in
    the background when that is smaller, so it takes its final value at construction.
    """
    cap: int = 2048
    bg_fasta: str | None = None
    bg_npz: str | None = None
    bg: BackgroundKmers | None = field(default=None, init=False, repr=False, compare=False)

    def __post_init__(self):
        if self.cap < 1 or self.cap > 2048:
            raise ValueError(f"cap must be between 1 and 2048, got {self.cap}")
        # frozen: the fields settled here are set through object.__setattr__, before
        # anything can have hashed the instance
        object.__setattr__(self, "bg", self._load_bg())
        if self.kmer != self.bg.kmer:
            raise ValueError(f"kmer mismatch: {self.kmer} != {self.bg.kmer}")
        max_counts = np.max(self.bg.counts)
        if max_counts < self.cap:
            object.__setattr__(self, "cap", int(max_counts))
        else:
            self.bg.counts = np.minimum(self.bg.counts, self.cap)

    def _load_bg(self) -> BackgroundKmers:
        """Load the background: from the FASTA if there is one, otherwise from the npz."""
        if self.bg_fasta is None and self.bg_npz is None:
            raise ValueError("No background provided: set either bg_fasta or bg_npz")
        if self.bg_fasta is not None:
            if self.bg_npz is not None:
                warnings.warn(
                    "both bg_fasta and bg_npz are set; bg_fasta takes priority and "
                    f"{self.bg_npz} is ignored",
                    stacklevel=4,
                )
            return BackgroundKmers.from_fasta(self.bg_fasta, self.kmer)
        return BackgroundKmers.load(self.bg_npz)

    def _kmer_scores(self, seq: str) -> tuple[np.ndarray, float]:
        return np.float16(self.bg.lookup(seq) / self.cap), float(self.cap)

@dataclass(frozen=True)
class SelfHomology(KmerHomology):

    def _kmer_scores(self, seq: str) -> tuple[np.ndarray, float]:
        repeats = BackgroundKmers.from_seqs([seq], self.kmer).lookup(seq) - 1
        multiplier = float(max(repeats.max(), 1))
        return np.float16(repeats / multiplier), multiplier

@dataclass(frozen=True)
class BoundaryMotif(AbstractScorer):
    """Overlaps that both open and close on a match of the motif score 0, the rest 1.

    The two ends are independent: the motif opening the overlap and the one closing it are
    separate occurrences. The motif is matched inside a lookahead so that overlapping
    occurrences are all found ("GG" starts at both 0 and 1 of "GGG"); note that only one
    match is taken per position, so a motif that can match several lengths at the same
    place ("[GC]{2,3}") contributes only its greedy end.
    """
    motif: str = "GC"

    def _fill(self, seq: str, grid: OverlapGrid, values: np.ndarray) -> float:
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

        # clamped only to keep the index in range: the cells whose overlap ends past the
        # sequence get a meaningless score here, which build_map overwrites with 1.
        shuts = closes[np.minimum(grid.ends, grid.seq_len)]

        bounded = opens[:, None] & shuts
        values[:] = (~bounded).astype(np.float16)
        return 1.0

