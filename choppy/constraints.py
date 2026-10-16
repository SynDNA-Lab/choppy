from dataclasses import dataclass, field
from abc import ABC, abstractmethod
import numpy as np
from Bio import SeqIO

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

@dataclass
class AbstractCondition(ABC):
    min_val: np.float16 | None = None
    max_val: np.float16 | None = None
    opt_val: np.float16 | None = None
    map: np.ndarray | None = field(default=None, init=False, repr=False)

    def build_map(self, seq: str, min_len: int, max_len: int) -> np.ndarray:
        self.map = np.full((len(seq) - min_len + 1, max_len - min_len + 1), np.nan, dtype=np.float16)
        self._fill_map(seq, min_len, max_len)
        return self.map

    @abstractmethod
    def _fill_map(self, seq: str, min_len: int, max_len: int) -> None: ...

class PolyX(AbstractCondition):

    def score(self, seq: str) -> np.float16:
        return np.float16(len(seq) > self.max_len)

    def _fill_map(self, seq: str, min_len: int, max_len: int) -> None:
        self.map[0, 0] = self.score(seq)

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

    def _fill_map(self, seq: str, min_len: int, max_len: int) -> None:
        scores = self._kmer_scores(seq)
        padded = np.concatenate([scores, np.full(max_len - min_len, np.nan, dtype=np.float16)])
        windows = np.lib.stride_tricks.sliding_window_view(padded, max_len - self.kmer + 1)
        self.map[:] = np.maximum.accumulate(windows, axis=1)[:, min_len - self.kmer:]

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
