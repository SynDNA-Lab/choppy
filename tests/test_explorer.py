"""Tests for choppy.explorer, run headless: values are typed into the fields and buttons
clicked as a user would, then the live objects and the figure are checked.

Skipped when the `gui` extra (ipywidgets, plotly, anywidget) is not installed.
"""
from dataclasses import dataclass

import numpy as np
import pytest

pytest.importorskip("ipywidgets")
pytest.importorskip("plotly")
pytest.importorskip("anywidget")

from choppy.constraints import AbstractScorer, BoundaryMotif, Condition, Relaxed  # noqa: E402
from choppy.explorer import (  # noqa: E402
    CRITICAL,
    MUTED,
    RELAXED,
    coverage,
    describe_relaxed,
    format_value,
    parse_value,
    relaxed_rows,
)
from choppy.fragmentor import Fragmentor  # noqa: E402

SEG = dict(min_segment_overlap_len=20, max_segment_overlap_len=30,
           min_segment_len=100, max_segment_len=200)
FRAG = dict(min_fragment_overlap_len=10, max_fragment_overlap_len=15,
            min_fragment_len=40, max_fragment_len=80)


def random_seq(n: int, seed: int = 0) -> str:
    return "".join(np.random.default_rng(seed).choice(list("ACGT"), n))


@dataclass(frozen=True)
class BlockRegion(AbstractScorer):
    """Scores 1 every overlap that touches [a, b), 0 the rest."""
    a: int
    b: int

    def _fill(self, seq, grid, values):
        values[:] = (grid.starts < self.b) & (grid.ends > self.a)
        return 0.0, 1.0


def fragmentor() -> Fragmentor:
    """1000 bp with three conditions: one blocking [300, 461), too long for any segment to
    cross; GC ends; and an objective."""
    f = Fragmentor(random_seq(1000), **SEG)
    f.add_condition(Condition(BlockRegion(300, 461), max_val=0.5, name="block"))
    f.add_condition(Condition(BoundaryMotif("[GC]"), max_val=0, name="gc"))
    f.add_condition(Condition(BoundaryMotif("A"), opt_val=0, name="ends in A"))
    return f


def traces(ex, block=None, kind=None) -> list:
    return [t for t in ex.fig.data
            if (block is None or t.meta["block"] == block) and (kind is None or t.meta["kind"] == kind)]


def lane(ex, block: int, y: int):
    """The stroke trace of one candidate lane of a block: pass 2, relaxed 1, manual 0."""
    found = [t for t in traces(ex, block, "lanes")
             if t.mode == "lines" and np.nanmax(np.asarray(t.y, dtype=float)) == y]
    return found[0] if found else None


def relaxed_lines(ex, block: int) -> list:
    return [t for t in traces(ex, block, "light") if t.name.startswith("relaxed")]


# --------------------------------------------------------------------------- helpers


@pytest.mark.parametrize(
    "text,kind,value",
    [("", float, None), (" none ", float, None), ("None", int, None), ("1.5", float, 1.5),
     ("1e3", int, 1000), ("-2", int, -2), ("inf", float, float("inf"))],
)
def test_parse_value(text, kind, value):
    assert parse_value(text, kind) == value


@pytest.mark.parametrize(
    "text,kind,match",
    [("abc", float, "not a number"), ("1.5", int, "whole number"), ("inf", int, "whole number")],
)
def test_parse_value_refuses_what_is_not_a_value(text, kind, match):
    with pytest.raises(ValueError, match=match):
        parse_value(text, kind)


def test_format_value():
    assert [format_value(v) for v in (None, np.int64(5), np.float16(0.5), 3.0, 1e-5)] == \
        ["", "5", "0.5", "3", "1e-05"]


def test_coverage_runs_are_the_covered_stretches():
    rng = np.random.default_rng(1)
    for _ in range(50):
        n = int(rng.integers(0, 30))
        starts = rng.integers(0, 200, n)
        ends = starts + rng.integers(1, 25, n)
        covered = np.zeros(260, dtype=bool)
        for s, e in zip(starts, ends):
            covered[s:e] = True  # touching intervals leave no hole, so they share a run
        edges = np.flatnonzero(np.diff(np.r_[0, covered.astype(int), 0]))
        expected = list(zip(edges[::2].tolist(), edges[1::2].tolist()))
        run_starts, run_ends, counts = coverage(starts, ends)
        assert list(zip(run_starts.tolist(), run_ends.tolist())) == expected
        assert counts.tolist() == [int(((starts >= a) & (starts < b)).sum()) for a, b in expected]


def test_relaxed_rows_are_its_area_within_the_conditions():
    f = Fragmentor(random_seq(1000), **SEG)
    c = Condition(BoundaryMotif("[GC]"), max_val=0, start=100, end=600)
    f.add_condition(c)
    c.add_relaxed(max_val=1, start=200, end=400)
    r = c.relaxed[0]
    assert relaxed_rows(c, r) == (200, 380)  # the last 20 bp overlap inside starts at 380
    c.start = 300  # the main area moved, and the relaxed version only applies inside it
    assert relaxed_rows(c, r) == (300, 380)
    c.start = 450
    assert relaxed_rows(c, r) is None


def test_describe_relaxed():
    assert describe_relaxed(Relaxed(0, None, 3.0, 2.0, 250, None)) == "#0 · max_val 3 · +2 · [250, …)"
    assert describe_relaxed(Relaxed(4, 0.5, None, 0.25, None, None)) == "#4 · min_val 0.5 · +0.25"


# --------------------------------------------------------------------------- opening


def test_explore_shows_the_first_four_conditions_by_default():
    f = Fragmentor(random_seq(1000), **SEG)
    for k in range(6):
        f.add_condition(Condition(BoundaryMotif("[GC]"), max_val=k, name=f"c{k}"))
    ex = f.explore()
    assert ex._blocks == [0, 1, 2, 3]
    assert [box.value for box in ex.checks] == [True] * 4 + [False] * 2
    assert len(traces(ex, kind="score")) == 4
    assert len(ex.panels.children) == 1 + 4  # the spacer, then one panel per block
    assert "application/vnd.jupyter.widget-view+json" in ex._repr_mimebundle_()


def test_explore_picks_the_shown_conditions_by_position_or_name():
    assert fragmentor().explore(show=["ends in A", 0])._blocks == [0, 2]
    with pytest.raises(ValueError, match="No condition at position 5"):
        fragmentor().explore(show=5)
    with pytest.raises(KeyError, match="nope"):
        fragmentor().explore(show="nope")


def test_explore_without_a_fragment_level_raises():
    with pytest.raises(ValueError, match="no fragment level"):
        fragmentor().explore(fragment=True)


def test_explore_the_fragment_level():
    f = Fragmentor(random_seq(1000), **SEG, **FRAG)
    f.add_condition(Condition(BoundaryMotif("[GC]"), max_val=0, name="gc"), fragment=True)
    ex = f.explore(fragment=True)
    assert ex.level is f.fragment_overlaps and (ex.min_len, ex.max_len) == (40, 80)
    assert "Fragment level" in ex.status.value


def test_explore_checks_its_view_and_row_height():
    f = fragmentor()
    with pytest.raises(ValueError, match="View range"):
        f.explore(a=50, b=50)
    with pytest.raises(ValueError, match="row_height"):
        f.explore(row_height=100)
    assert f.explore(a=100, b=300).fig.layout.xaxis.range == (99.5, 299.5)


def test_the_score_points_are_the_row_minima():
    f = fragmentor()
    (score,) = traces(f.explore(), block=1, kind="score")
    np.testing.assert_array_equal(score.y, f.get_condition("gc").row_minima())
    assert score.x[0] == 0 and score.x[-1] == f.segment_overlaps.grid.shape[0] - 1


def test_the_blocks_stack_at_fixed_heights():
    ex = fragmentor().explore(row_height=300)
    layout = ex.fig.layout
    assert layout.height == 3 * 300 + layout.margin.t + layout.margin.b
    tops = [layout[f"yaxis{2 * k + 1}" if k else "yaxis"].domain[1] for k in range(3)]
    # each score panel starts a title strip below its block's top, one block apart
    step = 300 / (3 * 300)
    assert tops == pytest.approx([1 - 24 / 900 - k * step for k in range(3)])
    assert all(p.layout.height == "300px" for p in ex.panels.children[1:])
    # a long panel or checklist scrolls inside its box (ipywidgets 8 has no overflow_y)
    assert all(p.layout.overflow == "hidden auto" for p in ex.panels.children[1:])
    assert ex.checklist.layout.overflow == "hidden auto"


# --------------------------------------------------------------------------- condition panels


def test_editing_a_threshold_updates_the_condition_and_its_line():
    f = fragmentor()
    ex = f.explore()
    ex.regenerate.click()
    c, panel = f.get_condition("gc"), ex.panel(1)
    panel.fields["max_val"].value = "1"
    assert c.max_val == 1.0
    assert [t.name for t in traces(ex, block=1, kind="light")] == ["max_val = 1"]
    assert f.segment_overlaps.candidates_stale and "stale" in ex.status.value
    assert lane(ex, 0, 2).opacity < 1  # the lanes fade until Regenerate


def test_a_refused_set_changes_nothing_until_the_other_field_makes_it_valid():
    f = fragmentor()
    ex = f.explore()
    c, panel = f.get_condition("gc"), ex.panel(1)
    panel.fields["min_val"].value = "0.5"  # above max_val 0
    assert c.min_val is None and c.version == 0
    assert "above max_val" in panel.message.value
    assert panel.fields["min_val"].value == "0.5"  # what was typed stays
    panel.fields["max_val"].value = "1"
    assert (c.min_val, c.max_val) == (0.5, 1.0) and panel.message.value == ""
    panel.fields["opt_val"].value = "abc"
    assert c.opt_val is None and "not a number" in panel.message.value


def test_an_area_past_the_sequence_is_cut_back_and_shown_so():
    f = fragmentor()
    ex = f.explore()
    panel = ex.panel(1)
    panel.fields["start"].value = "500"
    panel.fields["end"].value = "5000"
    c = f.get_condition("gc")
    assert (c.start, c.end) == (500, 1000)
    assert panel.fields["end"].value == "1000" and "past the sequence end" in panel.message.value
    washes = [(s.x0, s.x1) for s in ex.fig.layout.shapes if s.xref == "x3"]  # block 1's score axis
    assert washes == [(499.5, 980.5)]


def test_adding_and_removing_a_relaxed_version():
    f = fragmentor()
    ex = f.explore()
    c, panel = f.get_condition("block"), ex.panel(0)
    for name, value in (("max_val", "1"), ("penalty", "10"), ("start", "250"), ("end", "520")):
        panel.new[name].value = value
    panel.add.click()
    (r,) = c.relaxed
    assert (r.max_val, r.extra_penalty, r.start, r.end) == (1.0, 10.0, 250, 520)
    assert all(field.value == "" for field in panel.new.values())
    assert panel.relaxed.get_title(0) == "Relaxed (1)"
    (line,) = relaxed_lines(ex, 0)
    assert list(line.x) == [249.5, 500.5, None] and list(line.y) == [1.0, 1.0, None]
    assert list(line.text) == ["", f"#{r.id}", ""] and line.line.color == RELAXED[0]

    panel.relaxed_list.children[0].children[0].click()  # its remove button
    assert c.relaxed == () and relaxed_lines(ex, 0) == []
    assert panel.relaxed.get_title(0) == "Relaxed (0)"


def test_a_refused_relaxed_version_is_reported():
    f = fragmentor()
    ex = f.explore()
    panel = ex.panel(0)
    panel.new["max_val"].value = "0.25"  # lower than the main 0.5
    panel.add.click()
    assert f.get_condition("block").relaxed == ()
    assert "lower than main" in panel.message.value and panel.new["max_val"].value == "0.25"


def test_an_objective_cannot_be_relaxed():
    panel = fragmentor().explore().panel(2)
    assert panel.add.disabled and all(field.disabled for field in panel.new.values())
    assert "only a constraint" in panel.relaxed_note.value


def test_relaxed_versions_keep_their_colour_while_they_exist():
    f = fragmentor()
    ex = f.explore()
    c = f.get_condition("gc")
    ids = [c.add_relaxed(max_val=v) for v in (1, 2, 3, 4)]
    assert ex._relaxed_colors(1, c) == dict(zip(ids, [*RELAXED, MUTED]))  # never a fourth hue
    c.remove_relaxed(ids[0])
    assert ex._relaxed_colors(1, c) == dict(zip(ids[1:], [RELAXED[1], RELAXED[2], MUTED]))
    new = c.add_relaxed(max_val=5)
    assert ex._relaxed_colors(1, c)[new] == RELAXED[0]  # the slot freed by the removed one


# --------------------------------------------------------------------------- candidates


def test_regenerate_draws_the_lanes_and_the_gaps():
    f = fragmentor()
    ex = f.explore()
    assert "no candidates generated yet" in ex.status.value
    assert any("no candidates yet" in a.text for a in ex.fig.layout.annotations)
    ex.regenerate.click()
    ov = f.segment_overlaps
    assert ov.candidates is not None and not ov.candidates_stale

    gaps = ov.get_gapped_regions(200 + 1)
    assert gaps and ex._gaps == gaps and "gap" in ex.status.value
    washes = [(s.x0, s.x1) for s in ex.fig.layout.shapes if s.fillcolor == CRITICAL and s.xref == "x2"]
    assert washes == gaps

    c = ov.candidates
    run_starts, run_ends, _ = coverage(c.starts[1:-1], c.ends[1:-1])
    xs = np.asarray(lane(ex, 0, 2).x, dtype=float).reshape(-1, 4)  # start, middle, end, break
    np.testing.assert_array_equal(xs[:, 0], run_starts)
    np.testing.assert_array_equal(xs[:, 2], run_ends)
    assert lane(ex, 0, 1) is None and lane(ex, 0, 0) is None  # nothing relaxed, nothing by hand


def test_a_relaxation_fills_the_relaxed_lane_and_closes_the_gap():
    f = fragmentor()
    ex = f.explore()
    panel = ex.panel(0)
    for name, value in (("max_val", "1"), ("penalty", "10"), ("start", "250"), ("end", "520")):
        panel.new[name].value = value
    panel.add.click()
    ex.regenerate.click()
    assert ex._gaps == [] and "no gap is too long" in ex.status.value
    assert lane(ex, 0, 1) is not None


def test_adding_and_removing_a_candidate_by_hand():
    f = fragmentor()
    ex = f.explore()
    ex.regenerate.click()
    ex.new_start.value, ex.new_end.value, ex.new_penalty.value = "380", "405", "0.5"
    ex.add_candidate.click()
    ov = f.segment_overlaps
    assert list(ov.manual)[1:-1] == [(380, 405)] and (380, 405) in set(ov.candidates)
    assert ex.manual.get_title(0) == "Added by hand (1)"
    assert lane(ex, 0, 0) is not None and not ov.candidates_stale

    ex.manual_list.children[0].children[0].click()  # its remove button
    assert ov.manual.n_overlaps == 0 and (380, 405) not in set(ov.candidates)
    assert lane(ex, 0, 0) is None and ex.manual.get_title(0) == "Added by hand (0)"


def test_a_candidate_of_an_unusual_length_is_added_with_a_warning():
    f = fragmentor()
    ex = f.explore()
    ex.new_start.value, ex.new_end.value = "100", "160"  # 60 bp; the overlaps are 20-30
    ex.add_candidate.click()
    assert f.segment_overlaps.manual.n_overlaps == 1
    assert "overlap lengths 20–30; added anyway" in ex.message.value


@pytest.mark.parametrize(
    "start,end,penalty,match",
    [("500", "400", "", "start &lt; end"), ("", "400", "", "needs a start and an end"),
     ("100", "125", "-1", "non-negative"), ("1.5", "125", "", "whole number")],
)
def test_a_refused_candidate_is_reported(start, end, penalty, match):
    f = fragmentor()
    ex = f.explore()
    ex.new_start.value, ex.new_end.value, ex.new_penalty.value = start, end, penalty
    ex.add_candidate.click()
    assert f.segment_overlaps.manual.n_overlaps == 0 and match in ex.message.value


def test_editing_min_step():
    f = fragmentor()
    ex = f.explore()
    ex.regenerate.click()
    ov = f.segment_overlaps
    ex.min_step.value = "5"
    assert ov.min_step == 5 and ov.candidates_stale and "stale" in ex.status.value
    ex.min_step.value = "0"
    assert ov.min_step == 5 and "must be positive" in ex.message.value
    ex.min_step.value = "1"  # thins nothing: stored, and shown, as none
    assert ov.min_step is None and ex.min_step.value == ""


# --------------------------------------------------------------------------- checklist and refresh


def test_ticking_a_condition_adds_its_block_and_keeps_the_zoom():
    ex = fragmentor().explore(show=[0])
    ex.fig.layout.xaxis.range = (100, 200)  # as a zoom in the browser would
    ex.checks[1].value = True
    assert ex._blocks == [0, 1] and len(traces(ex, kind="score")) == 2
    assert ex.fig.layout.xaxis.range == (100, 200)
    ex.checks[0].value = False
    assert ex._blocks == [1] and len(ex.panels.children) == 2


def test_all_and_none():
    ex = fragmentor().explore(show=[0])
    ex.show_all.click()
    assert ex._blocks == [0, 1, 2] and all(box.value for box in ex.checks)
    ex.show_none.click()
    assert ex._blocks == [] and ex.fig is None and ex.panels.children == ()


def test_refresh_reads_changes_made_outside():
    f = fragmentor()
    ex = f.explore(show=[1])
    f.get_condition("gc").max_val = 2
    f.add_condition(Condition(BoundaryMotif("T"), max_val=0, name="late"))
    f.segment_overlaps.min_step = 7
    ex.refresh()
    assert ex.panel(1).fields["max_val"].value == "2"
    assert [t.name for t in traces(ex, block=1, kind="light")] == ["max_val = 2"]
    assert len(ex.checks) == 4 and ex.checks[3].description == "[3] late" and not ex.checks[3].value
    assert ex.min_step.value == "7"
