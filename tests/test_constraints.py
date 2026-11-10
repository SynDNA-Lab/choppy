"""Tests for choppy.constraints: the overlap grid, condition maps and the conditions.

Tests marked xfail document known bugs (see the reason). xfail is strict, so once a bug
is fixed the test starts to "XPASS" and fails the run: remove the marker then.
"""
import re
from dataclasses import dataclass

import numpy as np
import pytest

from choppy.constraints import (
    AbstractScorer,
    BackgroundKmers,
    BGHomology,
    BoundaryMotif,
    Condition,
    KmerHomology,
    MapStore,
    OverlapGrid,
    OverlapMap,
    SelfHomology,
    check_seq,
    encode_kmers,
    invalid_bases,
    reverse_complement,
    valid_kmers,
)

BASES = "ACGT"


def random_seq(n: int, seed: int = 0) -> str:
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list(BASES), n))


def kmer_code(kmer: str) -> int:
    code = 0
    for b in kmer:
        code = code * 4 + BASES.index(b)
    return code


# --------------------------------------------------------------------------- helpers


def test_reverse_complement():
    assert reverse_complement("AACGTT") == "AACGTT"
    assert reverse_complement("AAAC") == "GTTT"
    assert reverse_complement("") == ""


@pytest.mark.parametrize("kmer", [2, 4, 7, 16])
def test_encode_kmers_matches_base4(kmer):
    seq = random_seq(60, seed=kmer)
    enc = encode_kmers(seq, kmer)
    assert enc.dtype == np.uint64
    assert len(enc) == len(seq) - kmer + 1
    expected = [kmer_code(seq[p:p + kmer]) for p in range(len(seq) - kmer + 1)]
    assert enc.tolist() == expected


def test_check_seq_rejects_non_acgt():
    assert check_seq("ACGTACGT") == "ACGTACGT"
    with pytest.raises(ValueError, match="outside ACGT"):
        check_seq("ACGTNACGT")
    with pytest.raises(ValueError, match="position 0"):
        check_seq("acgt")


def test_invalid_bases_positions():
    assert invalid_bases("ACGT").tolist() == []
    assert invalid_bases("ANCGnT").tolist() == [1, 4]


@pytest.mark.parametrize(
    "seq,kmer,expected",
    [
        ("ACGTA", 3, [True, True, True]),
        ("ACNTA", 3, [False, False, False]),
        ("ACGTNACGT", 3, [True, True, False, False, False, True, True]),
        ("AC", 3, []),
    ],
)
def test_valid_kmers(seq, kmer, expected):
    assert valid_kmers(seq, kmer).tolist() == expected



# --------------------------------------------------------------------------- OverlapGrid


@pytest.fixture
def grid():
    return OverlapGrid(seq_len=50, min_len=5, max_len=9)


def test_grid_shape(grid):
    assert grid.shape == (46, 5)


def test_grid_cell_overlap_roundtrip_scalar(grid):
    assert grid.cell(10, 17) == (10, 2)
    assert grid.overlap(10, 2) == (10, 17)


def test_grid_cell_overlap_roundtrip_arrays(grid):
    i, j = np.meshgrid(np.arange(grid.shape[0]), np.arange(grid.shape[1]), indexing="ij")
    start, end = grid.overlap(i, j)
    i2, j2 = grid.cell(start, end)
    np.testing.assert_array_equal(i2, i)
    np.testing.assert_array_equal(j2, j)


@pytest.mark.parametrize("start,end", [(0, 4), (0, 10), (3, 3)])
def test_grid_cell_rejects_bad_length(grid, start, end):
    with pytest.raises(ValueError):
        grid.cell(start, end)


@pytest.mark.parametrize("i,j", [(-1, 0), (46, 0), (0, -1), (0, 5)])
def test_grid_overlap_rejects_out_of_bounds(grid, i, j):
    with pytest.raises(ValueError):
        grid.overlap(i, j)


@pytest.mark.parametrize("start", [-1, 46, 100])
def test_grid_cell_rejects_out_of_bounds_start(grid, start):
    with pytest.raises(ValueError, match="out of bounds"):
        grid.cell(start, start + 5)


def test_grid_cell_allows_end_past_sequence(grid):
    """The grid deliberately holds cells whose overlap runs past the end of the sequence;
    edge_mask is what drops them."""
    assert grid.cell(45, 54) == (45, 4)


def test_grid_broadcast_axes(grid):
    assert grid.starts.shape == (grid.shape[0], 1)
    assert grid.lengths.shape == (1, grid.shape[1])
    assert grid.ends.shape == grid.shape
    assert grid.starts[:, 0].tolist() == list(range(46))
    assert grid.lengths[0].tolist() == [5, 6, 7, 8, 9]
    i, j = 7, 3
    assert grid.ends[i, j] == grid.overlap(i, j)[1]


@pytest.mark.parametrize("min_edge_len", [19, 20, 30])
def test_grid_edge_mask_matches_brute_force(grid, min_edge_len):
    mask = grid.edge_mask(min_edge_len)
    assert mask.shape == grid.shape
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            start, end = grid.overlap(i, j)
            first_ok = end >= min_edge_len                    # [0, end) is long enough
            last_ok = grid.seq_len - start >= min_edge_len    # [start, L) is long enough
            assert mask[i, j] == (first_ok and last_ok), (start, end)


def test_grid_edge_mask_drops_overlaps_past_end_when_min_edge_exceeds_twice_max_len(grid):
    # The invariant Fragmentor._check_lengths relies on: min_edge_len > 2 * max_len
    mask = grid.edge_mask(2 * grid.max_len + 1)
    assert not np.any(mask & (grid.ends > grid.seq_len))


# --------------------------------------------------------------------------- OverlapMap


def test_overlap_map_lookups(grid):
    values = np.arange(np.prod(grid.shape)).reshape(grid.shape)
    m = OverlapMap(grid, values)
    assert m.get_cell(3, 2) == values[3, 2]
    assert m.get_overlap(3, 3 + 7) == values[3, 2]


@pytest.mark.parametrize("min_len,max_len", [(5, 9), (5, 5), (7, 8), (9, 9)])
def test_grid_slices_map_each_sub_cell_to_the_same_overlap(grid, min_len, max_len):
    sub = OverlapGrid(grid.seq_len, min_len, max_len)
    rows, cols = grid.slices(sub)
    np.testing.assert_array_equal(grid.starts[rows], sub.starts)
    np.testing.assert_array_equal(grid.ends[rows, cols], sub.ends)


@pytest.mark.parametrize(
    "sub,match",
    [
        (OverlapGrid(49, 5, 9), "sequence of 49"),
        (OverlapGrid(50, 4, 9), "not within"),
        (OverlapGrid(50, 6, 10), "not within"),
    ],
)
def test_grid_slices_rejects_a_grid_that_does_not_fit(grid, sub, match):
    with pytest.raises(ValueError, match=match):
        grid.slices(sub)


def test_overlap_map_view_shares_memory(grid):
    m = OverlapMap(grid, np.arange(np.prod(grid.shape), dtype=np.float16).reshape(grid.shape), 3.0)
    sub = OverlapGrid(grid.seq_len, 6, 8)
    v = m.view(sub)
    assert v.grid == sub and v.multiplier == 3.0
    assert np.shares_memory(v.values, m.values)
    assert v.get_overlap(10, 17) == m.get_overlap(10, 17)


# --------------------------------------------------------------------------- scorer and condition


@dataclass(frozen=True, eq=False)
class ConstScorer(AbstractScorer):
    """Fills the map from a fixed array, to test the base classes in isolation.

    Compared by identity, as the array cannot be hashed: two instances are two scorers.
    """
    fill: object = None
    __eq__ = object.__eq__
    __hash__ = object.__hash__

    def _fill(self, seq, grid, values):
        if self.fill is not None:
            values[:] = self.fill
        return 1.0


def bound(scorer, seq, grid, **thresholds) -> Condition:
    """A condition on scorer, bound to the map it builds over grid."""
    c = Condition(scorer, **thresholds)
    c.bind(scorer.build_map(seq, grid))
    return c


def test_abstract_scorer_cannot_be_instantiated():
    with pytest.raises(TypeError):
        AbstractScorer()


def test_scorers_are_frozen_and_hashed_by_parameters():
    assert BoundaryMotif("GC") == BoundaryMotif("GC")
    assert hash(BoundaryMotif("GC")) == hash(BoundaryMotif("GC"))
    assert BoundaryMotif("GC") != BoundaryMotif("CG")
    with pytest.raises(AttributeError):
        BoundaryMotif("GC").motif = "CG"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"min_val": 0.5},
        {"max_val": 0.5},
        {"opt_val": 0.0},
        {"opt_val": 0.0, "min_val": 0.1, "max_val": 0.9},
    ],
)
def test_condition_accepts_any_single_threshold(grid, kwargs):
    c = Condition(ConstScorer(), **kwargs)
    assert c.map is None and c.grid is None
    m = ConstScorer().build_map("A" * grid.seq_len, grid)
    c.bind(m)
    assert c.map is m and c.grid == grid


def test_condition_rejects_no_threshold():
    """With no threshold at all a condition is neither a constraint nor an objective,
    which is almost certainly a mistake in the caller."""
    with pytest.raises(ValueError, match="neither a constraint nor an objective"):
        Condition(ConstScorer())


def test_conditions_compare_by_identity():
    scorer = BoundaryMotif("GC")
    assert Condition(scorer, max_val=0) != Condition(scorer, max_val=0)


def test_build_map_checks_seq_length(grid):
    with pytest.raises(ValueError, match="long but the grid"):
        ConstScorer().build_map("A" * (grid.seq_len - 1), grid)


def test_build_map_defaults_to_worst_score(grid):
    m = ConstScorer().build_map("A" * grid.seq_len, grid)
    assert m.values.dtype == np.float16
    assert m.values.shape == grid.shape
    assert m.grid == grid
    assert np.all(m.values == 1.0)


def test_bind_rejects_a_second_map(grid):
    c = bound(ConstScorer(), "A" * grid.seq_len, grid, max_val=0.5)
    with pytest.raises(ValueError, match="already bound"):
        c.bind(ConstScorer().build_map("A" * grid.seq_len, grid))


def test_build_map_scores_overlaps_past_the_end_worst(grid):
    values = ConstScorer(fill=0.0).build_map("A" * grid.seq_len, grid).values
    past_end = grid.ends > grid.seq_len
    assert past_end.any() and (~past_end).any()
    assert np.all(values[past_end] == 1.0)
    assert np.all(values[~past_end] == 0.0)


def test_feasible_requires_map():
    with pytest.raises(ValueError, match="bound to a map"):
        Condition(ConstScorer(), max_val=0.5).feasible()


VALUES = np.array([[0.0, 0.25, 0.5], [0.75, 1.0, np.nan]], dtype=np.float16)
SMALL_GRID = OverlapGrid(seq_len=3, min_len=2, max_len=4)  # shape (2, 3)


@pytest.mark.parametrize(
    "kwargs,expected",
    [
        ({"max_val": 0.5}, [[1, 1, 1], [0, 0, 0]]),
        ({"min_val": 0.5}, [[0, 0, 1], [1, 1, 0]]),
        ({"min_val": 0.25, "max_val": 0.75}, [[0, 1, 1], [1, 0, 0]]),
        ({"opt_val": 0.0}, [[1, 1, 1], [1, 1, 1]]),
    ],
)
def test_feasible_thresholds(kwargs, expected):
    c = bound(ConstScorer(), "AAA", SMALL_GRID, **kwargs)
    c.map.values[:] = VALUES  # bypass the past-end rule: this tests the thresholds alone
    np.testing.assert_array_equal(c.feasible(), np.array(expected, dtype=bool))


# --------------------------------------------------------------------------- BackgroundKmers


@pytest.mark.parametrize("kmer", [1, 17])
def test_background_kmer_range(kmer):
    with pytest.raises(ValueError, match="between 2 and 16"):
        BackgroundKmers.from_seqs(["ACGTACGT"], kmer)


def test_background_counts_both_strands():
    bg = BackgroundKmers.from_seqs(["AAAAA"], 3)
    # 3 x AAA on the forward strand, 3 x TTT on the reverse strand
    assert bg.lookup("AAA").tolist() == [3]
    assert bg.lookup("TTT").tolist() == [3]
    assert bg.lookup("ACG").tolist() == [0]
    assert np.all(np.diff(bg.keys.astype(np.int64)) > 0)  # sorted, unique


def test_background_lookup_against_counter():
    seqs = [random_seq(200, seed=1), random_seq(150, seed=2)]
    k = 4
    bg = BackgroundKmers.from_seqs(seqs, k)
    allseqs = seqs + [reverse_complement(s) for s in seqs]
    query = random_seq(80, seed=3)
    expected = [
        sum(s[p:p + k] == query[q:q + k] for s in allseqs for p in range(len(s) - k + 1))
        for q in range(len(query) - k + 1)
    ]
    assert bg.lookup(query).tolist() == expected


def test_background_ignores_kmers_with_unknown_bases():
    with pytest.warns(UserWarning, match="outside ACGT"):
        bg = BackgroundKmers.from_seqs(["AAANAAA"], 3)
    # AAA at 0 and at 4 only; the three kmers covering the N are dropped on each strand
    assert bg.lookup("AAA").tolist() == [2]
    assert bg.lookup("TTT").tolist() == [2]
    assert len(bg.keys) == 2  # AAA and TTT


def test_background_uppercases_sequences():
    ref = BackgroundKmers.from_seqs(["ACGTACGT"], 4)
    low = BackgroundKmers.from_seqs(["acgtacgt"], 4)
    np.testing.assert_array_equal(low.keys, ref.keys)
    np.testing.assert_array_equal(low.counts, ref.counts)


def test_background_warns_and_skips_short_sequences():
    with pytest.warns(UserWarning, match="shorter than kmer"):
        bg = BackgroundKmers.from_seqs(["ACGTACGT", "AC"], 4)
    np.testing.assert_array_equal(bg.keys, BackgroundKmers.from_seqs(["ACGTACGT"], 4).keys)


def test_background_all_sequences_too_short():
    with pytest.raises(ValueError, match="as long as kmer"), pytest.warns(UserWarning):
        BackgroundKmers.from_seqs(["AC", "GT"], 4)


def test_background_no_valid_kmer():
    with pytest.raises(ValueError, match="No valid 3-mer"), pytest.warns(UserWarning):
        BackgroundKmers.from_seqs(["NNNNN"], 3)


def test_bg_homology_scores_unknown_background_region_as_zero(tmp_path):
    fasta = write_fasta(tmp_path / "bg.fa", ["AAAAAAAANNNNNNNN"])
    with pytest.warns(UserWarning, match="outside ACGT"):
        c = BGHomology(kmer=4, bg_fasta=fasta)
    assert c.bg.lookup("AAAA")[0] > 0
    assert c.bg.lookup("NNNN")[0] == 0  # the N kmer encodes as TTTT, which was dropped


def test_background_save_load_roundtrip(tmp_path):
    bg = BackgroundKmers.from_seqs([random_seq(100)], 5)
    path = tmp_path / "bg.npz"
    bg.save(str(path))
    loaded = BackgroundKmers.load(str(path))
    assert int(loaded.kmer) == 5
    np.testing.assert_array_equal(loaded.keys, bg.keys)
    np.testing.assert_array_equal(loaded.counts, bg.counts)


def test_background_from_fasta_uppercases(tmp_path):
    fasta = tmp_path / "bg.fa"
    fasta.write_text(">a\nacgtac\n>b\nGGGG\n")
    bg = BackgroundKmers.from_fasta(str(fasta), 3)
    ref = BackgroundKmers.from_seqs(["ACGTAC", "GGGG"], 3)
    np.testing.assert_array_equal(bg.keys, ref.keys)
    np.testing.assert_array_equal(bg.counts, ref.counts)


# --------------------------------------------------------------------------- KmerHomology


@dataclass(frozen=True, eq=False)
class StubHomology(KmerHomology):
    """KmerHomology with fixed per-kmer scores, to test the sliding-window max.

    Compared by identity, as the scores array cannot be hashed."""
    scores: object = None
    __eq__ = object.__eq__
    __hash__ = object.__hash__

    def _kmer_scores(self, seq):
        return np.asarray(self.scores, dtype=np.float16), 1.0


def brute_force_kmer_map(scores, grid, kmer):
    out = np.empty(grid.shape, dtype=np.float16)
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            s, e = grid.overlap(i, j)
            if e > grid.seq_len:
                out[i, j] = 1.0                      # runs past the end of the sequence
            elif e - s < kmer:
                out[i, j] = 0.0                      # holds no complete kmer: nothing against it
            else:
                out[i, j] = scores[s:e - kmer + 1].max()
    return out


@pytest.mark.parametrize("kmer,min_len,max_len", [(4, 4, 4), (4, 5, 9), (3, 8, 12), (5, 5, 20)])
def test_kmer_homology_is_window_max(kmer, min_len, max_len):
    L = 40
    rng = np.random.default_rng(kmer * 100 + min_len)
    scores = rng.random(L - kmer + 1).astype(np.float16) * 0.9
    grid = OverlapGrid(L, min_len, max_len)
    m = StubHomology(kmer=kmer, scores=scores).build_map("A" * L, grid)
    np.testing.assert_array_equal(m.values, brute_force_kmer_map(scores, grid, kmer))


@pytest.mark.parametrize("kmer,min_len,max_len", [(6, 5, 8), (9, 5, 12), (12, 4, 10)])
def test_kmer_homology_short_overlaps_pass(kmer, min_len, max_len):
    """An overlap shorter than kmer holds no complete kmer, so nothing counts against it."""
    L = 40
    rng = np.random.default_rng(kmer)
    scores = (rng.random(L - kmer + 1) * 0.9 + 0.1).astype(np.float16)  # all non-zero
    grid = OverlapGrid(L, min_len, max_len)
    m = StubHomology(kmer=kmer, scores=scores).build_map("A" * L, grid)
    short = m.values[:, :kmer - min_len]
    assert np.all(short[(grid.ends <= L)[:, :kmer - min_len]] == 0.0)
    np.testing.assert_array_equal(m.values, brute_force_kmer_map(scores, grid, kmer))


def test_kmer_homology_all_overlaps_shorter_than_kmer_pass():
    grid = OverlapGrid(40, 5, 8)
    m = StubHomology(kmer=20, scores=np.ones(21)).build_map("A" * 40, grid)
    inside = grid.ends <= grid.seq_len
    assert np.all(m.values[inside] == 0.0)
    assert np.all(m.values[~inside] == 1.0)
    assert m.multiplier == 1.0


def write_fasta(path, seqs):
    path.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(seqs)))
    return str(path)


def test_bg_homology_scores_are_capped_counts(tmp_path):
    fasta = write_fasta(tmp_path / "bg.fa", ["A" * 30])  # AAAA x27, TTTT x27
    c = BGHomology(kmer=4, cap=10, bg_fasta=fasta)
    assert c.cap == 10
    assert c.bg.counts.max() == 10
    scores, multiplier = c._kmer_scores("AAAAC")
    assert scores.dtype == np.float16
    assert scores.tolist() == [1.0, 0.0]
    assert multiplier == 10


def test_bg_homology_cap_lowered_to_max_count(tmp_path):
    fasta = write_fasta(tmp_path / "bg.fa", ["ACGTTGCA"])
    c = BGHomology(kmer=3, bg_fasta=fasta)
    assert c.cap == c.bg.counts.max()


def test_bg_homology_requires_background():
    with pytest.raises(ValueError, match="No background"):
        BGHomology(kmer=4)


def test_bg_homology_rejects_large_cap(tmp_path):
    fasta = write_fasta(tmp_path / "bg.fa", ["ACGTACGT"])
    with pytest.raises(ValueError, match="cap"):
        BGHomology(kmer=3, cap=4096, bg_fasta=fasta)


def test_bg_homology_kmer_mismatch_with_npz(tmp_path):
    npz = tmp_path / "bg.npz"
    BackgroundKmers.from_seqs(["ACGTACGTAC"], 5).save(str(npz))
    with pytest.raises(ValueError, match="kmer mismatch"):
        BGHomology(kmer=4, bg_npz=str(npz))


def test_bg_homology_from_npz_builds_map(tmp_path):
    npz = tmp_path / "bg.npz"
    BackgroundKmers.from_seqs([random_seq(300, seed=5)], 4).save(str(npz))
    seq = random_seq(60, seed=6)
    c = BGHomology(kmer=4, bg_npz=str(npz))
    grid = OverlapGrid(len(seq), 6, 10)
    m = c.build_map(seq, grid)
    scores, multiplier = c._kmer_scores(seq)
    expected = brute_force_kmer_map(scores, grid, 4)
    np.testing.assert_array_equal(m.values, expected)
    assert m.multiplier == multiplier == c.cap


def test_bg_homology_equality_ignores_the_loaded_background(tmp_path):
    """Equal source, kmer and cap make equal scorers, though each loaded its own copy."""
    fasta = write_fasta(tmp_path / "bg.fa", ["A" * 30])
    a, b = BGHomology(kmer=4, cap=10, bg_fasta=fasta), BGHomology(kmer=4, cap=10, bg_fasta=fasta)
    assert a.bg is not b.bg
    assert a == b and hash(a) == hash(b)
    assert a != BGHomology(kmer=4, cap=5, bg_fasta=fasta)


def test_bg_homology_lowered_cap_takes_part_in_equality(tmp_path):
    """cap is settled at construction, so caps that both clamp to the max count are equal."""
    fasta = write_fasta(tmp_path / "bg.fa", ["ACGTTGCA"])
    assert BGHomology(kmer=3, cap=100, bg_fasta=fasta) == BGHomology(kmer=3, cap=2048, bg_fasta=fasta)


def test_bg_homology_fasta_wins_over_npz(tmp_path):
    fasta = write_fasta(tmp_path / "bg.fa", ["AAAAAAAA"])
    npz = tmp_path / "bg.npz"
    BackgroundKmers.from_seqs(["ACGTACGTAC"], 4).save(str(npz))
    with pytest.warns(UserWarning, match="bg_fasta takes priority"):
        c = BGHomology(kmer=4, bg_fasta=fasta, bg_npz=str(npz))
    assert c.bg.lookup("AAAA")[0] > 0
    assert c.bg.lookup("ACGT")[0] == 0


@pytest.mark.parametrize("cap", [0, -1, 4096])
def test_bg_homology_rejects_bad_cap(tmp_path, cap):
    fasta = write_fasta(tmp_path / "bg.fa", ["ACGTACGT"])
    with pytest.raises(ValueError, match="cap must be between 1 and 2048"):
        BGHomology(kmer=3, cap=cap, bg_fasta=fasta)


# --------------------------------------------------------------------------- SelfHomology


def test_self_homology_unique_sequence_scores_zero():
    seq = "GTATCGGCTACCGCAAAAATAGTACCCTAT"  # every 5-mer unique on both strands
    assert BackgroundKmers.from_seqs([seq], 5).counts.max() == 1
    scores, multiplier = SelfHomology(kmer=5)._kmer_scores(seq)
    assert np.all(scores == 0)
    assert multiplier == 1


def test_self_homology_flags_repeats_and_inverted_repeats():
    unit = "GATTACA"
    seq = unit + "CCCGC" + reverse_complement(unit) + "GCGCC"
    scores, _ = SelfHomology(kmer=7)._kmer_scores(seq)
    assert scores[0] == 1.0                     # the repeat unit, via its reverse complement
    assert scores[len(unit) + 5] == 1.0         # and the reverse complement itself
    assert scores[1] == 0.0


def test_self_homology_penalises_palindromic_kmer():
    """Deliberate: a self-complementary kmer can anneal to itself and disturb the assembly,
    so it is counted on both strands and penalised even as a single occurrence."""
    seq = "ACGT" + "CAGGTCCTTGAGATTC"  # ACGT is its own reverse complement, occurs once
    scores, _ = SelfHomology(kmer=4)._kmer_scores(seq)
    assert scores[0] > 0.0


# --------------------------------------------------------------------------- BoundaryMotif


def brute_force_motif_map(seq, grid, motif):
    pattern = re.compile(motif)
    opens = {p for p in range(len(seq)) if pattern.match(seq, p)}
    closes = {m.end() for p in range(len(seq)) for m in [pattern.match(seq, p)] if m}
    out = np.ones(grid.shape, dtype=np.float16)
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            s, e = grid.overlap(i, j)
            if e <= grid.seq_len and s in opens and e in closes:
                out[i, j] = 0.0
    return out


def test_boundary_motif_non_overlapping_occurrences():
    seq = "AAGCAAAAGCAAAAAAAGCTTGCAAA"
    grid = OverlapGrid(len(seq), 4, 12)
    c = bound(BoundaryMotif(motif="GC"), seq, grid, max_val=0)
    expected = brute_force_motif_map(seq, grid, "GC")
    np.testing.assert_array_equal(c.map.values, expected)
    assert c.map.get_overlap(2, 10) == 0.0   # GCAAAAGC
    assert c.feasible().sum() == (expected == 0).sum()


def test_boundary_motif_regex():
    seq = "AAGCAAAACGAAAA"
    grid = OverlapGrid(len(seq), 4, 10)
    m = BoundaryMotif(motif="[GC][GC]").build_map(seq, grid)
    assert m.get_overlap(2, 10) == 0.0   # GC ... CG
    assert m.get_overlap(2, 9) == 1.0


def test_boundary_motif_overlapping_occurrences():
    seq = "AGGGAAAAAAGGGA"
    grid = OverlapGrid(len(seq), 4, 12)
    m = BoundaryMotif(motif="GG").build_map(seq, grid)
    expected = brute_force_motif_map(seq, grid, "GG")
    np.testing.assert_array_equal(m.values, expected)


def test_boundary_motif_past_end_scores_worst_when_sequence_ends_on_motif():
    seq = "ACGTACGTAA"
    grid = OverlapGrid(len(seq), 2, 4)
    m = BoundaryMotif(motif="A").build_map(seq, grid)
    np.testing.assert_array_equal(m.values, brute_force_motif_map(seq, grid, "A"))
    assert list(m.values[8]) == [0.0, 1.0, 1.0]   # only (8, 10) is a real overlap


# --------------------------------------------------------------------------- plotting


@pytest.fixture
def plt():
    """pyplot on the non-interactive backend; skips when the `plot` extra is not installed."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    yield plt
    plt.close("all")


def plotted_points(ax) -> tuple[np.ndarray, np.ndarray]:
    """x and y of the scatter plot_rows draws."""
    offsets = np.asarray(ax.collections[0].get_offsets())
    return offsets[:, 0], offsets[:, 1]


def test_plot_rows_takes_each_row_min_over_lengths(plt):
    grid = OverlapGrid(seq_len=10, min_len=2, max_len=4)  # shape (9, 3)
    fill = np.array(
        [[0.5, 0.25, 0.75], [1.0, 0.5, 0.125], [0.0, 1.0, 1.0], [0.75, 0.75, 0.75]] + [[1.0] * 3] * 5,
        dtype=np.float16,
    )
    c = bound(ConstScorer(fill=fill), "A" * 10, grid, max_val=0.5)
    x, y = plotted_points(c.plot_rows(0, 4))
    assert list(x) == [0, 1, 2, 3]
    assert list(y) == [0.25, 0.125, 0.0, 0.75]


@pytest.mark.parametrize(
    "kwargs,labels",
    [
        ({"max_val": 0.5}, ["max_val"]),
        ({"opt_val": 0.0}, ["opt_val"]),
        ({"min_val": 0.1, "max_val": 0.9, "opt_val": 0.2}, ["min_val", "max_val", "opt_val"]),
    ],
)
def test_plot_rows_draws_only_set_thresholds(plt, grid, kwargs, labels):
    c = bound(ConstScorer(fill=0.3), "A" * grid.seq_len, grid, **kwargs)
    ax = c.plot_rows(0, 2)
    drawn = [line.get_label().split(" =")[0] for line in ax.lines]
    assert drawn == labels
    for line, name in zip(ax.lines, labels):
        assert list(line.get_ydata()) == pytest.approx([kwargs[name]] * 2)


def test_plot_rows_uses_given_axes(plt, grid):
    c = bound(ConstScorer(), "A" * grid.seq_len, grid, max_val=0.5)
    _, ax = plt.subplots()
    assert c.plot_rows(0, 1, ax=ax) is ax


def test_plot_rows_titles_with_the_scorer(plt, grid):
    c = bound(BoundaryMotif(), "A" * grid.seq_len, grid, max_val=0.5)
    assert c.plot_rows(0, 1).get_title().startswith("BoundaryMotif:")


def test_plot_rows_titles_with_the_name_when_set(plt, grid):
    scorer = BoundaryMotif()
    c = Condition(scorer, max_val=0.5, name="gc ends")
    c.bind(scorer.build_map("A" * grid.seq_len, grid))
    assert c.plot_rows(0, 1).get_title().startswith("gc ends:")


def test_plot_rows_requires_map(plt):
    with pytest.raises(ValueError, match="bound to a map"):
        Condition(ConstScorer(), max_val=0.5).plot_rows(0, 1)


@pytest.mark.parametrize("a,b", [(-1, 2), (2, 2), (3, 1), (0, 10**6)])
def test_plot_rows_rejects_bad_range(plt, grid, a, b):
    c = bound(ConstScorer(), "A" * grid.seq_len, grid, max_val=0.5)
    with pytest.raises(ValueError, match="Row range"):
        c.plot_rows(a, b)


# --------------------------------------------------------------------------- shared maps


def view_scorers(tmp_path, seq):
    """One scorer of each kind, with a kmer that falls inside the union length range."""
    L = len(seq)
    fasta = write_fasta(tmp_path / "bg.fa", [random_seq(400, seed=11), seq[10:60]])
    scores = np.random.default_rng(3).random(L - 7 + 1).astype(np.float16)
    return [
        BoundaryMotif("[GC]"),
        StubHomology(kmer=7, scores=scores),
        SelfHomology(kmer=7),
        BGHomology(kmer=7, bg_fasta=fasta),
    ]


@pytest.mark.parametrize("min_len,max_len", [(4, 6), (8, 12), (15, 20), (4, 20)])
def test_view_of_a_union_map_equals_a_standalone_map(tmp_path, min_len, max_len):
    """Every score depends on the overlap alone, not on the grid it is laid out on, so a
    level reading a slice of a wider map sees exactly the map it would build itself."""
    seq = random_seq(80, seed=12)
    union = OverlapGrid(len(seq), 4, 20)
    sub = OverlapGrid(len(seq), min_len, max_len)
    for scorer in view_scorers(tmp_path, seq):
        own = scorer.build_map(seq, sub)
        view = scorer.build_map(seq, union).view(sub)
        np.testing.assert_array_equal(view.values, own.values, err_msg=type(scorer).__name__)
        assert view.multiplier == own.multiplier


def test_map_store_builds_each_scorer_once(grid):
    store = MapStore("A" * grid.seq_len, grid)
    scorer = BoundaryMotif("A")
    assert scorer not in store
    m = store.get(scorer)
    assert scorer in store and len(store) == 1
    assert store.get(scorer) is m
    assert store.get(BoundaryMotif("A")) is m        # equal, not identical
    assert store.get(BoundaryMotif("AA")) is not m
    assert len(store) == 2


def test_map_store_views_share_one_map(grid):
    store = MapStore(random_seq(grid.seq_len), grid)
    scorer = BoundaryMotif("[GC]")
    a = store.view(scorer, OverlapGrid(grid.seq_len, 5, 7))
    b = store.view(BoundaryMotif("[GC]"), OverlapGrid(grid.seq_len, 6, 9))
    assert len(store) == 1
    assert np.shares_memory(a.values, b.values)


def test_map_store_maps_are_read_only(grid):
    store = MapStore("A" * grid.seq_len, grid)
    m = store.get(BoundaryMotif("A"))
    with pytest.raises(ValueError, match="read-only"):
        m.values[0, 0] = 0.5
    with pytest.raises(ValueError, match="read-only"):
        store.view(BoundaryMotif("A"), OverlapGrid(grid.seq_len, 6, 7)).values[0, 0] = 0.5


def test_map_store_checks_seq_length(grid):
    with pytest.raises(ValueError, match="long but the grid"):
        MapStore("A" * (grid.seq_len + 1), grid)
