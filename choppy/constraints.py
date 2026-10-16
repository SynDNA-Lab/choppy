from dataclasses import dataclass, field
from abc import ABC, abstractmethod
import numpy as np
from Bio import SeqIO

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

    @abstractmethod
    def score(self, seq: str) -> np.float16: ...

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
        lut = np.full(256, 255, np.uint8)
        for i, b in enumerate(b"ACGT"): lut[b] = i
        all_codes = []
        for seq in sequences:
            codes = lut[np.frombuffer(seq.encode(), np.uint8)]
            m = len(codes) - kmer + 1
            enc = np.zeros(m, dtype=np.uint64)
            for j in range(kmer):
                enc = (enc << np.uint64(2)) | codes[j:j + m]
            all_codes.append(enc)

        keys, counts = np.unique(np.concatenate(all_codes), return_counts=True)
        return keys, counts

@dataclass
class bg_homology(AbstractCondition):
    bg_fasta: str | None = None
    kmer: int = 15
    trie: dict | None = field(default=None, init=False, repr=False)
    multiplier: float = field(init=False, repr=False)

    def score(self, seq: str) -> np.float16:
        return np.float16(len(seq) > self.max_len)

    def _fill_map(self, seq: str, min_len: int, max_len: int) -> None:
        bg = BackgroundKmers.from_fasta(self.bg_fasta, self.kmer)
        self.map[0, 0] = self.score(seq)