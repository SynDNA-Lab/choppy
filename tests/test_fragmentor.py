import itertools
import numpy as np
import pytest

from choppy.constraints import AbstractCondition, BoundaryMotif, OverlapGrid
from choppy.fragmentor import Fragmentor, Overlaps


class ConstCondition(AbstractCondition):
    def __init__(self, fill=None, **kwargs):
        super().__init__(**kwargs)
        self.fill = fill

    def _fill_map(self, seq):
        if self.fill is not None:
            self.map.values[:] = self.fill


def random_seq(n: int, seed: int = 0) -> str:
    return "".join(np.random.default_rng(seed).choice(list("ACGT"), n))


# --------------------------------------------------------------------------- Overlaps


def test_overlaps_grid_and_shape():
    ov = Overlaps("A" * 30, min_len=4, max_len=7)
    assert ov.grid == OverlapGrid(30, 4, 7)
    assert ov.shape == (27, 4)


def test_overlaps_min_step_one_disables_thinning():
    assert Overlaps("A" * 30, 4, 7, min_step=1).min_step is None


def test_overlaps_instances_do_not_share_lists():
    a, b = Overlaps("A" * 30, 4, 7), Overlaps("A" * 30, 4, 7)
    a.add_condition(ConstCondition(max_val=0.5))
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
    c = ConstCondition(**kwargs)
    ov.add_condition(c)
    assert (c in ov.objectives) == is_objective
    assert (c in ov.constraints) == is_constraint
    assert c.grid == ov.grid
    assert c.map is not None and c.map.values.shape == ov.shape


def test_add_condition_rejects_meaningless_condition():
    """A condition with no threshold is neither objective nor constraint: add_condition
    builds its map first, so the rejection has to leave both lists untouched."""
    ov = Overlaps("A" * 30, 4, 7)
    c = ConstCondition()
    with pytest.raises(ValueError, match="neither a constraint nor an objective"):
        ov.add_condition(c)
    assert ov.objectives == [] and ov.constraints == []
    assert c.map is None


def test_fragmentor_add_condition_rejects_meaningless_condition():
    f = Fragmentor("ACGT" * 100, **SEG, **FRAG)
    for fragment in (False, True):
        with pytest.raises(ValueError, match="neither a constraint nor an objective"):
            f.add_condition(BoundaryMotif(), fragment=fragment)
    assert f.segment_overlaps.constraints == [] and f.fragment_overlaps.constraints == []


def test_feasible_mask_without_constraints_is_all_true():
    ov = Overlaps("A" * 30, 4, 7)
    assert ov.feasible_mask().all()


def test_feasible_mask_is_and_of_constraints():
    ov = Overlaps("A" * 30, 4, 7)
    rng = np.random.default_rng(0)
    a = rng.random(ov.shape).astype(np.float16)
    b = rng.random(ov.shape).astype(np.float16)
    ov.add_condition(ConstCondition(fill=a, max_val=0.6))
    ov.add_condition(ConstCondition(fill=b, min_val=0.3))
    ov.add_condition(ConstCondition(fill=np.ones(ov.shape), opt_val=0.0))  # objective only
    np.testing.assert_array_equal(ov.feasible_mask(), (a <= 0.6) & (b >= 0.3))


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
    ov.add_condition(BoundaryMotif(motif="[GC]", max_val=0))
    cands = ov.generate_candidates(min_edge_len=40)
    assert cands == ov.candidates
    assert cands == sorted(cands)
    assert cands, "expected some candidates with a single-base motif"
    for s, e in cands:
        assert 10 <= e - s <= 15
        assert e >= 40 and len(seq) - s >= 40
        assert seq[s] in "GC" and seq[e - 1] in "GC"
    expected = ov.feasible_mask() & ov.grid.edge_mask(40)
    assert len(cands) == expected.sum()


def test_generate_candidates_empty():
    ov = Overlaps("A" * 100, 10, 15)
    ov.add_condition(BoundaryMotif(motif="GC", max_val=0))
    assert ov.generate_candidates(min_edge_len=31) == []


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
    seg_c, frag_c = BoundaryMotif(max_val=0), BoundaryMotif(max_val=0)
    f.add_condition(seg_c)
    f.add_condition(frag_c, fragment=True)
    assert f.segment_overlaps.constraints == [seg_c]
    assert f.fragment_overlaps.constraints == [frag_c]
    assert seg_c.grid.min_len == 20 and frag_c.grid.min_len == 10


def test_fragmentor_same_condition_on_both_levels_is_rejected():
    f = Fragmentor("ACGT" * 100, **SEG, **FRAG)
    c = BoundaryMotif(max_val=0)
    f.add_condition(c)
    with pytest.raises(ValueError, match="already bound"):
        f.add_condition(c, fragment=True)


def test_fragmentor_fragment_condition_without_fragment_level():
    f = Fragmentor("ACGT" * 100, **SEG)
    with pytest.raises(ValueError, match="no fragment level"):
        f.add_condition(BoundaryMotif(max_val=0), fragment=True)


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
