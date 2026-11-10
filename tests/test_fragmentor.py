import itertools
from dataclasses import dataclass

import numpy as np
import pytest

from choppy.constraints import AbstractScorer, BoundaryMotif, Condition, MapStore, OverlapGrid
from choppy.fragmentor import Candidates, Fragmentor, GapError, NoAssemblyError, Overlaps


@dataclass(frozen=True, eq=False)
class ConstScorer(AbstractScorer):
    """Fills the map from a fixed array; compared by identity, as arrays cannot be hashed."""
    fill: object = None
    __eq__ = object.__eq__
    __hash__ = object.__hash__

    def _fill(self, seq, grid, values):
        if self.fill is not None:
            values[:] = self.fill
        return 1.0


def const(fill=None, **thresholds) -> Condition:
    return Condition(ConstScorer(fill), **thresholds)


def motif(motif="GC", **thresholds) -> Condition:
    return Condition(BoundaryMotif(motif), **thresholds)


def random_seq(n: int, seed: int = 0) -> str:
    return "".join(np.random.default_rng(seed).choice(list("ACGT"), n))


# --------------------------------------------------------------------------- Overlaps


def test_overlaps_grid_and_shape():
    ov = Overlaps("A" * 30, min_len=4, max_len=7)
    assert ov.grid == OverlapGrid(30, 4, 7)
    assert ov.shape == (27, 4)


def test_overlaps_min_step_one_disables_thinning():
    assert Overlaps("A" * 30, 4, 7, min_step=1).min_step is None


def test_overlaps_without_a_store_gets_its_own():
    ov = Overlaps("A" * 30, 4, 7)
    assert ov.maps.grid == ov.grid and len(ov.maps) == 0
    ov.add_condition(motif("A", max_val=0))
    assert len(ov.maps) == 1


def test_overlaps_reads_a_view_of_a_wider_store():
    seq = random_seq(60)
    store = MapStore(seq, OverlapGrid(60, 3, 12))
    ov = Overlaps(seq, 5, 8, maps=store)
    c = motif("[GC]", max_val=0)
    ov.add_condition(c)
    assert c.grid == ov.grid
    assert np.shares_memory(c.map.values, store.get(BoundaryMotif("[GC]")).values)
    np.testing.assert_array_equal(c.map.values, BoundaryMotif("[GC]").build_map(seq, ov.grid).values)


@pytest.mark.parametrize(
    "store,match",
    [
        (MapStore("A" * 30, OverlapGrid(30, 4, 6)), "not within"),
        (MapStore("C" * 30, OverlapGrid(30, 4, 7)), "different sequence"),
    ],
)
def test_overlaps_rejects_a_store_it_does_not_fit(store, match):
    with pytest.raises(ValueError, match=match):
        Overlaps("A" * 30, 4, 7, maps=store)


def test_overlaps_instances_do_not_share_lists():
    a, b = Overlaps("A" * 30, 4, 7), Overlaps("A" * 30, 4, 7)
    a.add_condition(const(max_val=0.5))
    assert b.constraints == [] and b.objectives == []


@pytest.mark.parametrize(
    "kwargs,is_objective,is_constraint",
    [
        ({"opt_val": 0.0}, True, False),
        ({"max_val": 0.5}, False, True),
        ({"min_val": 0.5}, False, True),
        ({"opt_val": 0.0, "max_val": 0.5}, True, True),
    ],
)
def test_add_condition_routing(kwargs, is_objective, is_constraint):
    ov = Overlaps("A" * 30, 4, 7)
    c = const(**kwargs)
    ov.add_condition(c)
    assert (c in ov.objectives) == is_objective
    assert (c in ov.constraints) == is_constraint
    assert c.grid == ov.grid
    assert c.map is not None and c.map.values.shape == ov.shape


def test_feasible_mask_without_constraints_is_all_true():
    ov = Overlaps("A" * 30, 4, 7)
    assert ov.feasible_mask().all()


def test_feasible_mask_is_and_of_constraints():
    ov = Overlaps("A" * 30, 4, 7)
    rng = np.random.default_rng(0)
    a = rng.random(ov.shape).astype(np.float16)
    b = rng.random(ov.shape).astype(np.float16)
    ov.add_condition(const(fill=a, max_val=0.6))
    ov.add_condition(const(fill=b, min_val=0.3))
    ov.add_condition(const(fill=np.ones(ov.shape), opt_val=0.0))  # objective only
    inside = ov.grid.ends <= ov.grid.seq_len  # past the end every score is 1
    expected = np.where(inside, (a <= 0.6) & (b >= 0.3), (1.0 <= 0.6) & (1.0 >= 0.3))
    np.testing.assert_array_equal(ov.feasible_mask(), expected)


# ------------------------------------------------------------------ thinning


def kept_overlaps(ov, mask):
    s, e = ov.grid.overlap(*np.nonzero(mask))
    return list(zip(s.tolist(), e.tolist()))


def manhattan(a, b):
    return abs(a[0] - b[0]) + abs(a[1] - b[1])


def test_thin_without_min_step_is_identity():
    ov = Overlaps("A" * 40, 5, 12)
    mask = np.random.default_rng(1).random(ov.shape) < 0.5
    np.testing.assert_array_equal(ov._thin(mask.copy()), mask)


@pytest.mark.parametrize("step", [2, 3, 5])
def test_thin_returns_subset_and_keeps_first(step):
    ov = Overlaps("A" * 40, 5, 12, min_step=step)
    mask = np.random.default_rng(step).random(ov.shape) < 0.6
    out = ov._thin(mask.copy())
    assert not np.any(out & ~mask)
    first = tuple(np.argwhere(mask)[0])
    assert out[first]


def test_thin_suppresses_diagonal_neighbour():
    # (0, 7) and (1, 7): start moves by 1, end stays -> distance 1 < 2
    ov = Overlaps("A" * 40, 5, 12, min_step=2)
    mask = np.zeros(ov.shape, bool)
    mask[ov.grid.cell(0, 7)] = True
    mask[ov.grid.cell(1, 7)] = True
    assert kept_overlaps(ov, ov._thin(mask)) == [(0, 7)]


def test_thin_keeps_far_neighbour():
    # (0, 7) and (1, 8): distance 2 == min_step -> both kept
    ov = Overlaps("A" * 40, 5, 12, min_step=2)
    mask = np.zeros(ov.shape, bool)
    mask[ov.grid.cell(0, 7)] = True
    mask[ov.grid.cell(1, 8)] = True
    assert kept_overlaps(ov, ov._thin(mask)) == [(0, 7), (1, 8)]


def test_thin_does_not_wrap_negative_columns():
    ov = Overlaps("A" * 40, 5, 8, min_step=2)            # 4 length columns
    mask = np.zeros(ov.shape, bool)
    mask[ov.grid.cell(0, 5)] = True                      # row 0, col 0
    mask[ov.grid.cell(1, 9)] = True                      # row 1, last col: (1,9) vs (0,5) is 5 apart
    assert kept_overlaps(ov, ov._thin(mask)) == [(0, 5), (1, 9)]


def test_thin_suppresses_same_start_neighbour():
    ov = Overlaps("A" * 40, 5, 12, min_step=2)
    mask = np.zeros(ov.shape, bool)
    mask[ov.grid.cell(0, 7)] = True
    mask[ov.grid.cell(0, 8)] = True
    assert kept_overlaps(ov, ov._thin(mask)) == [(0, 7)]


@pytest.mark.parametrize("step", [2, 3, 4, 7])
def test_thin_result_respects_min_step(step):
    """Documented contract: kept overlaps are pairwise >= min_step apart (start + end)."""
    ov = Overlaps("A" * 40, 5, 12, min_step=step)
    mask = np.random.default_rng(step).random(ov.shape) < 0.6
    kept = kept_overlaps(ov, ov._thin(mask.copy()))
    too_close = [(a, b) for a, b in itertools.combinations(kept, 2) if manhattan(a, b) < step]
    assert too_close == []
    # and nothing was dropped without a kept overlap close to it
    for cand in kept_overlaps(ov, mask):
        assert cand in kept or any(manhattan(cand, k) < step for k in kept)


# ------------------------------------------------------------------ candidates


def test_generate_candidates_respects_constraints_and_edges():
    seq = random_seq(200, seed=3)
    ov = Overlaps(seq, 10, 15)
    ov.add_condition(motif("[GC]", max_val=0))
    cands = ov.generate_candidates(min_edge_len=40)
    assert cands is ov.candidates
    assert list(cands) == sorted(cands)
    assert np.all(cands.extra_penalty == 0)
    pairs = list(cands)
    assert pairs[0] == (0, 0) and pairs[-1] == (len(seq), len(seq))
    assert cands.n_overlaps > 0, "expected some candidates with a single-base motif"
    for s, e in pairs[1:-1]:
        assert 10 <= e - s <= 15
        assert e >= 40 and len(seq) - s >= 40
        assert seq[s] in "GC" and seq[e - 1] in "GC"
    expected = ov.feasible_mask() & ov.grid.edge_mask(40)
    assert cands.n_overlaps == expected.sum()


def test_generate_candidates_empty():
    ov = Overlaps("A" * 100, 10, 15)
    ov.add_condition(motif("GC", max_val=0))
    cands = ov.generate_candidates(min_edge_len=31)
    assert cands.n_overlaps == 0 and list(cands) == [(0, 0), (100, 100)]


def test_generate_candidates_thinned():
    seq = random_seq(200, seed=4)
    plain = Overlaps(seq, 10, 15).generate_candidates(40)
    thinned = Overlaps(seq, 10, 15, min_step=4).generate_candidates(40)
    assert set(thinned) < set(plain)


# --------------------------------------------------------------------------- Fragmentor

SEG = dict(min_segment_overlap_len=20, max_segment_overlap_len=30,
           min_segment_len=100, max_segment_len=200)
FRAG = dict(min_fragment_overlap_len=10, max_fragment_overlap_len=15,
            min_fragment_len=40, max_fragment_len=80)


def test_fragmentor_segment_only():
    f = Fragmentor("acgt" * 100, **SEG)
    assert f.seq == "ACGT" * 100
    assert f.seq_len == 400
    assert f.use_fragments is False
    assert f.fragment_overlaps is None
    assert f.segment_overlaps.grid == OverlapGrid(400, 20, 30)
    assert f.segment_overlaps.seq == f.seq
    assert f.grid == OverlapGrid(400, 20, 30)
    assert f.segment_overlaps.maps is f.maps


def test_fragmentor_grid_spans_both_levels():
    f = Fragmentor("ACGT" * 100, **SEG, **FRAG)
    assert f.grid == OverlapGrid(400, 10, 30)
    assert f.maps.grid == f.grid
    assert f.segment_overlaps.maps is f.maps and f.fragment_overlaps.maps is f.maps


def test_fragmentor_grid_follows_the_lowered_max_overlap():
    with pytest.warns(UserWarning, match="Lowering it to 29"):
        f = Fragmentor("ACGT" * 100, **{**SEG, "min_segment_len": 60})
    assert f.grid == OverlapGrid(400, 20, 29)


def test_fragmentor_levels_share_one_map_per_scorer():
    seq = random_seq(400, seed=7)
    f = Fragmentor(seq, **SEG, **FRAG)
    seg_c, frag_c = motif("[GC]", max_val=0), motif("[GC]", min_val=1)
    f.add_condition(seg_c)
    f.add_condition(frag_c, fragment=True)
    assert len(f.maps) == 1
    shared = f.maps.get(BoundaryMotif("[GC]")).values
    for c, overlaps in ((seg_c, f.segment_overlaps), (frag_c, f.fragment_overlaps)):
        # the two levels' columns are disjoint here, so each view shares only with the full map
        assert np.shares_memory(c.map.values, shared)
        own = BoundaryMotif("[GC]").build_map(seq, overlaps.grid)
        np.testing.assert_array_equal(c.map.values, own.values)
    f.add_condition(motif("GC", max_val=0))
    assert len(f.maps) == 2


def test_fragmentor_warns_when_most_shared_lengths_are_unused():
    # segment overlaps 40-42 and fragment overlaps 10-11: 5 lengths used out of 10-42
    seg = {**SEG, "min_segment_overlap_len": 40, "max_segment_overlap_len": 42}
    frag = {**FRAG, "max_fragment_overlap_len": 11}
    with pytest.warns(UserWarning, match="28 of those 33 lengths are used by neither level"):
        Fragmentor("ACGT" * 100, **seg, **frag)


def test_fragmentor_with_fragments():
    f = Fragmentor("ACGT" * 100, **SEG, **FRAG, min_step=3)
    assert f.use_fragments is True
    assert f.fragment_overlaps.grid == OverlapGrid(400, 10, 15)
    assert f.fragment_overlaps.min_step == 3
    assert f.segment_overlaps.min_step == 3


@pytest.mark.parametrize("missing", list(FRAG))
def test_fragmentor_partial_fragment_params(missing):
    partial = {k: v for k, v in FRAG.items() if k != missing}
    with pytest.raises(ValueError, match="all set or all None"):
        Fragmentor("ACGT" * 100, **SEG, **partial)


@pytest.mark.parametrize(
    "override,match",
    [
        ({"min_segment_overlap_len": 2.0}, "must be int"),
        ({"min_segment_overlap_len": True}, "must be int"),
        ({"min_segment_overlap_len": 0}, "must be positive"),
        ({"min_segment_overlap_len": 31}, r"min_segment_overlap_len \(31\) > max"),
        ({"min_segment_len": 61, "max_segment_len": 61}, "max_segment_len"),
        ({"max_segment_len": 99}, "max_segment_len"),
    ],
)
def test_fragmentor_segment_length_validation(override, match):
    with pytest.raises(ValueError, match=match):
        Fragmentor("ACGT" * 100, **{**SEG, **override})


def test_fragmentor_minimal_valid_lengths():
    Fragmentor("ACGT" * 100, **{**SEG, "min_segment_len": 61, "max_segment_len": 62})


@pytest.mark.parametrize(
    "override,match",
    [
        ({"min_fragment_overlap_len": -1}, "must be positive"),
        ({"max_fragment_len": 40}, "max_fragment_len"),
    ],
)
def test_fragmentor_fragment_length_validation(override, match):
    with pytest.raises(ValueError, match=match):
        Fragmentor("ACGT" * 100, **SEG, **{**FRAG, **override})


def test_fragmentor_add_condition_routing():
    f = Fragmentor("ACGT" * 100, **SEG, **FRAG)
    seg_c, frag_c = motif(max_val=0), motif(max_val=0)
    f.add_condition(seg_c)
    f.add_condition(frag_c, fragment=True)
    assert f.segment_overlaps.constraints == [seg_c]
    assert f.fragment_overlaps.constraints == [frag_c]
    assert seg_c.grid.min_len == 20 and frag_c.grid.min_len == 10


def test_fragmentor_same_condition_on_both_levels_is_rejected():
    f = Fragmentor("ACGT" * 100, **SEG, **FRAG)
    c = motif(max_val=0)
    f.add_condition(c)
    with pytest.raises(ValueError, match="already bound"):
        f.add_condition(c, fragment=True)


def test_fragmentor_fragment_condition_without_fragment_level():
    f = Fragmentor("ACGT" * 100, **SEG)
    with pytest.raises(ValueError, match="no fragment level"):
        f.add_condition(motif(max_val=0), fragment=True)


@pytest.mark.parametrize("n", [1, 16, 99])
def test_fragmentor_rejects_too_short_sequence(n):
    with pytest.raises(ValueError, match="shorter than min_segment_len"):
        Fragmentor("A" * n, **SEG)


def test_fragmentor_accepts_sequence_exactly_min_segment_len():
    f = Fragmentor("A" * SEG["min_segment_len"], **SEG)
    assert f.seq_len == SEG["min_segment_len"]


def test_fragmentor_rejects_non_acgt_sequence():
    with pytest.raises(ValueError, match="outside ACGT"):
        Fragmentor("ACGN" * 100, **SEG)


@pytest.mark.parametrize("step", [0, -3])
def test_fragmentor_rejects_non_positive_min_step(step):
    with pytest.raises(ValueError, match="min_step must be positive"):
        Fragmentor("ACGT" * 100, **SEG, min_step=step)


def test_fragmentor_min_step_one_disables_thinning():
    f = Fragmentor("ACGT" * 100, **SEG, min_step=1)
    assert f.segment_overlaps.min_step is None


@pytest.mark.parametrize(
    "override,match",
    [
        ({"min_fragment_len": 150, "max_fragment_len": 180},
         r"min_fragment_len \(150\) must be <= min_segment_len"),
    ],
)
def test_fragmentor_fragment_min_len_must_fit_inside_segment(override, match):
    with pytest.raises(ValueError, match=match):
        Fragmentor("ACGT" * 100, **SEG, **{**FRAG, **override})


# ------------------------------------------------- the smaller element's max is clamped


def overlap_bounds(overlaps):
    return overlaps.grid.min_len, overlaps.grid.max_len


@pytest.mark.parametrize(
    "min_segment_len,expected_max_overlap",
    [
        (60, 29),   # 2 * 30 == 60, so the longest overlap is 29
        (61, 30),   # exactly at the limit: nothing to lower
        (41, 20),   # lowered right down to min_segment_overlap_len
    ],
)
def test_segment_overlap_max_clamped_to_fit_segment(min_segment_len, expected_max_overlap):
    params = {**SEG, "min_segment_len": min_segment_len}
    if expected_max_overlap == SEG["max_segment_overlap_len"]:
        f = Fragmentor("ACGT" * 100, **params)
    else:
        with pytest.warns(UserWarning, match="Lowering it to"):
            f = Fragmentor("ACGT" * 100, **params)
    assert overlap_bounds(f.segment_overlaps) == (SEG["min_segment_overlap_len"], expected_max_overlap)
    assert min_segment_len > 2 * f.segment_overlaps.grid.max_len


def test_fragment_overlap_max_clamped_to_fit_fragment():
    with pytest.warns(UserWarning, match=r"max_fragment_overlap_len \(15\) is too long"):
        f = Fragmentor("ACGT" * 100, **SEG, **{**FRAG, "min_fragment_len": 30})
    assert overlap_bounds(f.fragment_overlaps) == (10, 14)
    assert overlap_bounds(f.segment_overlaps) == (20, 30)  # the other level is untouched


def test_max_fragment_len_clamped_to_max_segment_len():
    with pytest.warns(UserWarning, match=r"max_fragment_len \(250\) is longer than"):
        f = Fragmentor("ACGT" * 100, **SEG, **{**FRAG, "max_fragment_len": 250})
    assert f.max_fragment_len == SEG["max_segment_len"] == 200
    assert f.max_fragment_len > f.min_fragment_len  # the clamp cannot invert the range


def test_clamp_keeps_fragment_range_valid_at_the_extreme():
    """min_fragment_len == min_segment_len is the tightest the clamp can get squeezed."""
    tight = {**FRAG, "min_fragment_len": SEG["min_segment_len"], "max_fragment_len": 9999}
    with pytest.warns(UserWarning, match="Lowering it to"):
        f = Fragmentor("ACGT" * 100, **SEG, **tight)
    assert f.min_fragment_len == 100 and f.max_fragment_len == 200


@pytest.mark.parametrize(
    "params,match",
    [
        ({**SEG, "min_segment_len": 41, "min_segment_overlap_len": 21},
         r"min_segment_overlap_len \(21\) leaves no room"),
        ({**SEG, "min_segment_len": 31, "min_segment_overlap_len": 16},
         r"leaves no room in a min_segment_len \(31\) segment"),
    ],
)
def test_min_overlap_len_with_no_room_is_an_error(params, match):
    """Nothing can be lowered to rescue a minimal overlap that does not fit."""
    with pytest.raises(ValueError, match=match):
        Fragmentor("ACGT" * 100, **params)


def test_no_warning_for_consistent_lengths(recwarn):
    Fragmentor("ACGT" * 100, **SEG, **FRAG)
    assert [w for w in recwarn.list if issubclass(w.category, UserWarning)] == []


def test_fragmentor_validates_once_for_every_grid():
    """Fragmentor is the only validator, so its checks must leave both grids usable:
    max_len <= seq_len on each level, which is what the condition maps assume."""
    f = Fragmentor("ACGT" * 100, **SEG, **FRAG)
    for overlaps in (f.segment_overlaps, f.fragment_overlaps):
        assert overlaps.grid.max_len <= overlaps.grid.seq_len
        assert overlaps.grid.shape[0] > 0 and overlaps.grid.shape[1] > 0


# --------------------------------------------------------------------------- gap check


@dataclass(frozen=True)
class BlockRegion(AbstractScorer):
    """Scores 1 every overlap that touches [a, b), 0 the rest."""
    a: int
    b: int

    def _fill(self, seq, grid, values):
        values[:] = (grid.starts < self.b) & (grid.ends > self.a)
        return 1.0


def block(a, b) -> Condition:
    """Rejects every overlap that touches [a, b)."""
    return Condition(BlockRegion(a, b), max_val=0.5)


@pytest.mark.parametrize("a", [300, 500])
def test_chop_reports_unbridgeable_gap(a):
    # SEG: the shortest crossing opens on a 20 bp overlap ending at a and closes on a
    # 20 bp one starting at a + 161, so it is 161 + 2 * 20 = 201 bp, one over 200
    f = Fragmentor(random_seq(1000), **SEG)
    f.add_condition(block(a, a + 161))
    with pytest.raises(GapError, match="add candidate overlaps manually") as info:
        f.chop()
    assert info.value.gaps == [(a - 20, a + 181)]
    assert f"[{a - 20}, {a + 181}): 201 bp" in str(info.value)
    assert isinstance(info.value, NoAssemblyError)


def test_chop_bridges_gap_just_under_the_limit():
    f = Fragmentor(random_seq(1000), **SEG)
    f.add_condition(block(500, 660))
    asm = f.chop()
    assert asm.segments[0][0] == 0 and asm.segments[-1][1] == 1000


def test_chop_reports_whole_sequence_without_candidates():
    f = Fragmentor(random_seq(1000), **SEG)
    f.add_condition(const(fill=1.0, max_val=0.5))
    with pytest.raises(GapError) as info:
        f.chop()
    assert info.value.gaps == [(0, 1000)]


def test_gap_lengths_include_the_edges():
    ov = Overlaps("A" * 100, 5, 8)
    ov.candidates = Candidates(100, [30], [36])
    assert ov.gap_lengths().tolist() == [36, 70]
    assert ov.get_gapped_regions(1) == [(0, 36), (30, 100)]
    assert ov.get_gapped_regions(50) == [(30, 100)]


def test_gap_lengths_with_unsorted_ends():
    # B starts before C but ends after it, so the split between A and B is crossed
    # most cheaply by A -> C (75 bp), not A -> B (90 bp)
    ov = Overlaps("A" * 1000, 15, 40)
    ov.candidates = Candidates(1000, [100, 150, 155], [130, 190, 175])
    assert ov.gap_lengths().tolist() == [130, 75, 25, 845]
    assert ov.get_gapped_regions(76) == [(0, 130), (155, 1000)]


def test_gap_check_allows_a_gap_of_exactly_max_len():
    # every crossing is exactly max_segment_len (200), which a segment may be
    f = Fragmentor("A" * 1000, **SEG)
    f.segment_overlaps.candidates = Candidates(1000, [80, 260, 440, 620, 800],
                                               [100, 280, 460, 640, 820])
    f.gap_check()
    f.segment_overlaps.candidates = Candidates(1000, [80, 260, 441, 620, 800],
                                               [100, 280, 461, 640, 820])
    with pytest.raises(GapError) as info:
        f.gap_check()
    assert info.value.gaps == [(260, 461)]


def test_gap_lengths_without_candidates():
    ov = Overlaps("A" * 100, 5, 8)
    ov.candidates = Candidates(100)
    assert ov.gap_lengths().tolist() == [100]
    assert ov.get_gapped_regions(1) == [(0, 100)]
    assert ov.get_gapped_regions(101) == []


def test_chop_short_sequence_needs_no_candidates():
    f = Fragmentor(random_seq(150), **SEG)
    f.add_condition(const(fill=1.0, max_val=0.5))
    assert f.chop().segments == [(0, 150)]


# --------------------------------------------------------------------------- Candidates


def test_candidates_sorted_by_start_then_end():
    c = Candidates(100, [5, 1, 5, 3], [9, 4, 7, 8], [0.5, 0.1, 0.2, 0.3])
    assert list(c) == [(0, 0), (1, 4), (3, 8), (5, 7), (5, 9), (100, 100)]
    # moved with their overlaps
    assert c.extra_penalty.tolist() == [0.0, 0.1, 0.3, 0.2, 0.5, 0.0]


def test_candidates_default_and_scalar_penalty():
    assert Candidates(100, [1, 2], [5, 6]).extra_penalty.tolist() == [0.0] * 4
    assert Candidates(100, [1, 2], [5, 6], 2.0).extra_penalty.tolist() == [0.0, 2.0, 2.0, 0.0]


def test_candidates_empty():
    c = Candidates(100)
    assert c.n_overlaps == 0 and list(c) == [(0, 0), (100, 100)]


def test_candidates_add_keeps_order():
    c = Candidates(100, [10, 30], [20, 40])
    c.add(25, 35, extra_penalty=1.0)
    c.add([0, 10, 90], [15, 18, 100])
    assert list(c) == [(0, 0), (0, 15), (10, 18), (10, 20), (25, 35), (30, 40), (90, 100),
                       (100, 100)]
    assert c.extra_penalty.tolist() == [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]


def test_candidates_add_skips_duplicates_by_default():
    c = Candidates(100, [10], [20], 1.0)
    c.add([10, 12], [20, 22], 5.0)
    assert list(c)[1:-1] == [(10, 20), (12, 22)]
    assert c.extra_penalty.tolist()[1:-1] == [1.0, 5.0]


def test_candidates_add_updates_duplicates():
    c = Candidates(100, [10], [20], 1.0)
    c.add([10, 12], [20, 22], 5.0, skip_duplicates=False)
    assert list(c)[1:-1] == [(10, 20), (12, 22)]
    assert c.extra_penalty.tolist()[1:-1] == [5.0, 5.0]


@pytest.mark.parametrize("skip,expected", [(True, 1.0), (False, 3.0)])
def test_candidates_duplicates_within_a_batch(skip, expected):
    c = Candidates(100)
    c.add([7, 7, 7], [9, 9, 9], [1.0, 2.0, 3.0], skip_duplicates=skip)
    assert list(c)[1:-1] == [(7, 9)]
    assert c.extra_penalty.tolist()[1:-1] == [expected]


def test_candidates_arrays_are_read_only():
    c = Candidates(100, [1], [5])
    for arr in (c.starts, c.ends, c.extra_penalty):
        with pytest.raises(ValueError):
            arr[0] = 0


@pytest.mark.parametrize(
    "args,match",
    [
        (([1, 2], [5]), "equal length"),
        (([5], [5]), "start < end"),
        (([-1], [5]), "start < end"),
        (([90], [101]), "end <= seq_len"),
        (([0], [0]), "start < end"),  # sentinels cannot be added
        (([1, 2], [5, 6], [1.0, 2.0, 3.0]), "extra_penalty"),
    ],
)
def test_candidates_rejects_bad_input(args, match):
    with pytest.raises(ValueError, match=match):
        Candidates(100, *args)


def test_chop_keeps_existing_candidates():
    f = Fragmentor(random_seq(1000), **SEG)
    f.generate_candidates()
    f.segment_overlaps.candidates.add(500, 525, extra_penalty=1.0)
    f.chop()
    assert (500, 525) in set(f.segment_overlaps.candidates)


# --------------------------------------------------------------------------- names and bookkeeping


def test_constraints_and_objectives_follow_the_conditions():
    ov = Overlaps("A" * 30, 4, 7)
    both, only_opt = const(max_val=0.5, opt_val=0.0), const(opt_val=0.0)
    ov.add_condition(both)
    ov.add_condition(only_opt)
    assert ov.conditions == [both, only_opt]
    assert ov.constraints == [both]
    assert ov.objectives == [both, only_opt]


def test_get_condition_by_name():
    f = Fragmentor("ACGT" * 100, **SEG, **FRAG)
    gc, at = motif("GC", max_val=0, name="gc"), motif("AT", max_val=0, name="at")
    f.add_condition(gc)
    f.add_condition(at, fragment=True)
    assert f.get_condition("gc") is gc
    assert f.get_condition("at", fragment=True) is at
    with pytest.raises(KeyError, match=r"No condition named 'at'.*\['gc'\]"):
        f.get_condition("at")


def test_get_condition_without_fragment_level():
    with pytest.raises(ValueError, match="no fragment level"):
        Fragmentor("ACGT" * 100, **SEG).get_condition("gc", fragment=True)


def test_duplicate_name_on_a_level_is_rejected_before_binding():
    ov = Overlaps("A" * 30, 4, 7)
    ov.add_condition(motif("GC", max_val=0, name="x"))
    dup = motif("AT", max_val=0, name="x")
    with pytest.raises(ValueError, match="already a condition named 'x'"):
        ov.add_condition(dup)
    assert dup.map is None and len(ov.conditions) == 1


def test_same_name_on_both_levels_is_fine():
    f = Fragmentor("ACGT" * 100, **SEG, **FRAG)
    f.add_condition(motif(max_val=0, name="gc"))
    f.add_condition(motif(max_val=0, name="gc"), fragment=True)


def test_second_condition_on_an_equal_scorer_is_rejected_before_binding():
    ov = Overlaps("A" * 30, 4, 7)
    ov.add_condition(motif("GC", max_val=0, name="strict"))
    second = motif("GC", max_val=1)  # an equal scorer, though a separate instance
    with pytest.raises(ValueError, match=r"already a condition on BoundaryMotif.*\(strict\)"):
        ov.add_condition(second)
    assert second.map is None and len(ov.conditions) == 1
    ov.add_condition(motif("CG", max_val=1))  # a different scorer is fine


def test_label_is_the_name_or_the_scorer_class():
    assert motif(max_val=0).label == "BoundaryMotif"
    assert motif(max_val=0, name="gc").label == "gc"


@pytest.fixture
def plt():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    yield plt
    plt.close("all")


def plotted_overlaps():
    ov = Overlaps(random_seq(60), 5, 8)
    ov.add_condition(motif("GC", max_val=0, name="gc"))
    ov.add_condition(motif("AT", max_val=0))
    ov.add_condition(motif("CC", opt_val=0, name="cc"))  # objective only
    return ov


def test_overlaps_plot_rows_all_constraints_by_default(plt):
    axes = plotted_overlaps().plot_rows()
    assert [ax.get_title() for ax in axes] == ["[0] gc", "[1] BoundaryMotif"]


def test_overlaps_plot_rows_by_name_and_index(plt):
    axes = plotted_overlaps().plot_rows(["cc", 1])
    assert [ax.get_title() for ax in axes] == ["cc", "[1] BoundaryMotif"]


def test_overlaps_plot_rows_unknown_name(plt):
    with pytest.raises(KeyError, match="No condition named 'nope'"):
        plotted_overlaps().plot_rows("nope")


def test_fragmentor_plot_rows_by_name(plt):
    f = Fragmentor(random_seq(400), **SEG, **FRAG)
    f.add_condition(motif("GC", max_val=0, name="gc"), fragment=True)
    axes = f.plot_rows(fragment=True, indices="gc", a=0, b=50)
    assert [ax.get_title() for ax in axes] == ["gc"]
