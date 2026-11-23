"""Interactive explorer of one level's conditions and candidates: Fragmentor.explore.

A notebook widget, ipywidgets around one plotly FigureWidget, that edits the live objects.
Each shown condition is a block of the figure:
- its row minima, as Condition.plot_rows draws them;
- its thresholds and area;
- the bounds of its relaxed versions, over the starts they apply to;
- below these, lanes with the level's candidates and the gaps between them.

Next to each block is a panel where the condition's thresholds, area and relaxed versions
can be edited. A toolbar holds the condition checklist, min_step, candidate regeneration and
the overlaps added by hand.

Edits never regenerate the candidates: they show as stale until Regenerate is pressed
(chop and gap_check regenerate stale candidates anyway).

All the blocks are one figure on purpose. Plotly draws scattergl with WebGL, which takes
about two browser contexts per figure, and browsers keep only about 16, so a figure per
condition would start going blank after a handful of conditions. The panels have fixed
heights so that they stay level with their blocks.

Needs the `gui` extra: ipywidgets, plotly and anywidget (plotly's FigureWidget runs on it).
"""
from __future__ import annotations

import html
import time
import warnings
from typing import TYPE_CHECKING

import numpy as np

try:
    import anywidget  # noqa: F401  (plotly's FigureWidget runs on it)
    import ipywidgets as widgets
    import plotly.graph_objects as go
except ImportError as error:  # pragma: no cover - depends on what is installed
    raise ImportError(
        "The explorer needs the `gui` extra (ipywidgets, plotly and anywidget): "
        "`poetry install --all-extras`, or `pip install 'choppy[gui]'`"
    ) from error

from choppy.constraints import THRESHOLDS

if TYPE_CHECKING:
    from choppy.constraints import Condition, Relaxed
    from choppy.fragmentor import Overlaps

DEFAULT_SHOWN = 4     # conditions plotted at first, so that the first view is not overplotted
MIN_ROW_HEIGHT = 170  # px: a block's fixed strips plus a usable score panel
MAX_LISTED = 200      # hand-added candidates listed with a remove button

# The data hues are categorical slots validated together, all pairs, on the light surface.
# A hue is never generated or cycled: past three relaxed versions of one condition, the
# rest are MUTED, and each still carries its #id label.
SURFACE, GRID, AXIS = "#fcfcfb", "#e1e0d9", "#c3c2b7"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
SCORE, OPT, AREA = "#2a78d6", "#008300", "#eb6834"
RELAXED = ("#e87ba4", "#eda100", "#1baf7a")
CRITICAL = "#d03b3b"  # gaps and errors, always with a label
FONT = "system-ui, -apple-system, 'Segoe UI', sans-serif"
LANES = ("pass", "relaxed", "manual")  # the candidate lanes, top to bottom
MARGIN_LEFT, MARGIN_RIGHT = 72, 12
PANEL_WIDTH = 340


# --------------------------------------------------------------------------- helpers


def parse_value(text: str, kind: type = float):
    """The value typed into a field: None when it is empty or "none", otherwise a number.

    An int field takes whole numbers only ("1e3" is 1000). Raises ValueError otherwise.
    """
    s = text.strip()
    if s.lower() in ("", "none"):
        return None
    try:
        value = float(s)
    except ValueError:
        raise ValueError(f"{text!r} is not a number") from None
    if kind is int:
        if not value.is_integer():
            raise ValueError(f"{text!r} is not a whole number")
        return int(value)
    return value


def format_value(value) -> str:
    """A stored value as a field shows it, with "" for None."""
    if value is None:
        return ""
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return f"{float(value):g}"


def coverage(starts, ends) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The union of the intervals [start, end) as runs sorted by start: the run starts, the
    run ends, and the number of intervals in each run. Intervals that overlap or touch
    share a run."""
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    if starts.size == 0:
        return (np.zeros(0, dtype=np.int64),) * 3
    order = np.argsort(starts, kind="stable")
    s, e = starts[order], ends[order]
    reach = np.maximum.accumulate(e)
    first = np.flatnonzero(np.r_[True, s[1:] > reach[:-1]])
    last = np.r_[first[1:] - 1, len(s) - 1]
    return s[first], reach[last], last - first + 1


def relaxed_rows(condition: Condition, r: Relaxed) -> tuple[int, int] | None:
    """The overlap starts where relaxed version r of condition can change anything
    (Condition.relaxed_rows), or None when no overlap fits there."""
    return condition.relaxed_rows(r)


def describe_relaxed(r: Relaxed) -> str:
    """One line for a relaxed version: its id, the bounds it sets, its penalty and its area
    (… for an open side)."""
    parts = [f"#{r.id}"]
    for name, value in (("min_val", r.min_val), ("max_val", r.max_val)):
        if value is not None:
            parts.append(f"{name} {format_value(value)}")
    parts.append(f"+{r.extra_penalty:g}")
    if r.start is not None or r.end is not None:
        lo = "…" if r.start is None else f"{r.start:,}"
        hi = "…" if r.end is None else f"{r.end:,}"
        parts.append(f"[{lo}, {hi})")
    return " · ".join(parts)


def _attempt(action) -> tuple[str | None, list[str]]:
    """Run action, catching what a refused edit raises: (error or None, warnings)."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            action()
            error = None
        except (ValueError, TypeError, KeyError) as e:
            error = e.args[0] if isinstance(e, KeyError) and e.args else str(e)
    return error, [str(w.message) for w in caught]


def _notice(error: str | None = None, warned=(), info: str | None = None) -> str:
    """HTML for a message line: an error, warnings and/or a plain note, each marked."""
    lines = []
    if error:
        lines.append(f"<span style='color:{CRITICAL}'>✖ {html.escape(error)}</span>")
    lines += [f"<span style='color:{INK2}'>⚠ {html.escape(w)}</span>" for w in warned]
    if info:
        lines.append(f"<span style='color:{INK2}'>{html.escape(info)}</span>")
    return "<br>".join(lines)


def _text(description: str, placeholder: str = "none", width: int = 150, label: int = 52):
    return widgets.Text(description=description, placeholder=placeholder, continuous_update=False,
                        layout=widgets.Layout(width=f"{width}px"),
                        style={"description_width": f"{label}px"})


def _button(description: str, tooltip: str = "", width: str = "auto"):
    return widgets.Button(description=description, tooltip=tooltip,
                          layout=widgets.Layout(width=width))


# --------------------------------------------------------------------------- figure layout


class _Geometry:
    """Where each block goes. Every shown condition gets row_height px of the plot area,
    split from the top into a title strip, its score panel, a gap, its candidate lanes, and
    a space holding the lanes' x tick labels."""
    TITLE, GAP, LANES, SPACE = 24, 6, 46, 26
    TOP, BOTTOM = 6, 44  # figure margins; the bottom one holds the x axis title

    def __init__(self, n_blocks: int, row_height: int):
        self.n, self.h = n_blocks, row_height
        self.plot = n_blocks * row_height

    @property
    def height(self) -> int:
        return self.plot + self.TOP + self.BOTTOM

    def _y(self, px: float) -> float:
        """Paper y of the point px below the top of the plot area."""
        return 1 - px / self.plot

    def top(self, k: int) -> float:
        return self._y(k * self.h)

    def score(self, k: int) -> list[float]:
        bottom = (k + 1) * self.h - self.SPACE - self.LANES - self.GAP
        return [self._y(bottom), self._y(k * self.h + self.TITLE)]

    def lanes(self, k: int) -> list[float]:
        bottom = (k + 1) * self.h - self.SPACE
        return [self._y(bottom), self._y(bottom - self.LANES)]


def _axes(k: int) -> tuple[str, str, str, str]:
    """The trace axes of block k: score x, score y, lanes x, lanes y."""
    def name(axis: str, i: int) -> str:
        return axis if i == 1 else f"{axis}{i}"
    s, l = 2 * k + 1, 2 * k + 2
    return name("x", s), name("y", s), name("x", l), name("y", l)


def _key(axis: str) -> str:
    """The layout key of a trace axis: "x3" -> "xaxis3"."""
    return f"{axis[0]}axis{axis[1:]}"


def _legend(k: int) -> str:
    return "legend" if k == 0 else f"legend{k + 1}"


def _block_layout(geo: _Geometry, k: int, name: str) -> dict:
    """The axes and the legend of block k; its x axes follow the first block's."""
    sx, sy, lx, ly = _axes(k)
    ticks = dict(size=10, color=MUTED)
    x = dict(domain=[0, 1], showgrid=True, gridcolor=GRID, gridwidth=1, zeroline=False,
             showline=True, linecolor=AXIS, linewidth=1, ticks="", tickfont=ticks,
             tickformat=",d", automargin=False)
    last = k == geo.n - 1
    return {
        _key(sx): dict(x, anchor=sy, showticklabels=False, **({} if k == 0 else {"matches": "x"})),
        _key(sy): dict(domain=geo.score(k), anchor=sx, gridcolor=GRID, gridwidth=1, zeroline=False,
                       showline=True, linecolor=AXIS, ticks="", tickfont=ticks, automargin=False,
                       title=dict(text="score", font=dict(size=10, color=MUTED), standoff=6)),
        _key(lx): dict(x, anchor=ly, matches="x", showticklabels=True,
                       **({"title": dict(text=f"overlap start (bp) · {name} level",
                                         font=dict(size=11, color=INK2), standoff=4)} if last else {})),
        _key(ly): dict(domain=geo.lanes(k), anchor=lx, range=[-0.6, 2.6], fixedrange=True,
                       tickmode="array", tickvals=[2, 1, 0], ticktext=list(LANES), showgrid=False,
                       zeroline=False, showline=False, ticks="", automargin=False,
                       tickfont=dict(size=9, color=MUTED)),
        _legend(k): dict(x=1, xanchor="right", y=geo.top(k), yanchor="top", orientation="h",
                         bgcolor="rgba(252,252,251,0.85)", font=dict(size=10, color=INK2),
                         itemsizing="constant"),
    }


# --------------------------------------------------------------------------- the explorer


class _ConditionPanel:
    """The editable parameters of one condition, level with its block."""

    def __init__(self, explorer: LevelExplorer, index: int, condition: Condition, height: int):
        self.explorer, self.index, self.condition = explorer, index, condition
        self._syncing = False

        self.title = widgets.HTML()
        self.fields = {name: _text(name) for name in THRESHOLDS}
        for f in self.fields.values():
            f.observe(self._on_field, names="value")
        self.message = widgets.HTML()

        self.relaxed_list = widgets.VBox()
        self.relaxed_note = widgets.HTML()
        self.new = {
            "min_val": _text("min_val"), "max_val": _text("max_val"),
            "penalty": _text("penalty", placeholder="0"),
            "start": _text("start"), "end": _text("end"),
        }
        self.add = _button("Add relaxed", "add a relaxed version with these bounds, penalty and area")
        self.add.on_click(lambda _: explorer._add_relaxed(self))
        form = widgets.VBox([
            widgets.HBox([self.new["min_val"], self.new["max_val"]]),
            widgets.HBox([self.new["penalty"]]),
            widgets.HBox([self.new["start"], self.new["end"]]),
            self.add,
        ])
        self.relaxed = widgets.Accordion(children=[widgets.VBox([self.relaxed_list, self.relaxed_note, form])],
                                         selected_index=None)

        f = self.fields
        self.widget = widgets.VBox(
            [self.title,
             widgets.HBox([f["min_val"], f["max_val"]]),
             widgets.HBox([f["opt_val"]]),
             widgets.HBox([f["start"], f["end"]]),
             self.relaxed, self.message],
            layout=widgets.Layout(height=f"{height}px", overflow="hidden auto", margin="0",
                                  padding="0 0 0 10px", border_bottom=f"1px solid {GRID}"),
        )
        self.sync()

    def sync(self) -> None:
        """Show the condition as it is now."""
        c = self.condition
        self._syncing = True
        try:
            for name, field_ in self.fields.items():
                field_.value = format_value(getattr(c, name))
        finally:
            self._syncing = False
        kinds = [k for k, on in (("constraint", c.is_constraint), ("objective", c.is_objective)) if on]
        self.title.value = (
            f"<b>[{self.index}] {html.escape(c.label)}</b> "
            f"<span style='color:{MUTED}'>{' · '.join(kinds)}</span><br>"
            f"<span style='color:{MUTED};font-size:11px;word-break:break-all'>"
            f"{html.escape(repr(c.scorer))}</span>"
        )
        self.sync_relaxed()

    def sync_relaxed(self) -> None:
        c = self.condition
        colors = self.explorer._relaxed_colors(self.index, c)
        rows = []
        for r in sorted(c.relaxed, key=lambda r: r.id):
            note = ""
            if relaxed_rows(c, r) is None:
                note = " <i>(no effect: outside the condition's area)</i>"
            elif r.min_val is None and r.max_val is None:
                note = " <i>(no effect: no looser bound)</i>"
            remove = _button("✕", f"remove relaxed #{r.id}", width="32px")
            remove.on_click(lambda _, rid=r.id: self.explorer._remove_relaxed(self, rid))
            text = widgets.HTML(f"<span style='color:{colors[r.id]}'>━━</span> "
                                f"{html.escape(describe_relaxed(r))}{note}")
            rows.append(widgets.HBox([remove, text]))
        self.relaxed_list.children = rows or [widgets.HTML("<i>none</i>")]
        self.relaxed.set_title(0, f"Relaxed ({len(c.relaxed)})")
        for w in (*self.new.values(), self.add):
            w.disabled = not c.is_constraint
        self.relaxed_note.value = ("" if c.is_constraint else
                                   "<i>only a constraint (min_val and/or max_val) can be relaxed</i>")

    def say(self, error: str | None = None, warned=(), info: str | None = None) -> None:
        self.message.value = _notice(error, warned, info)

    def _on_field(self, change) -> None:
        if not self._syncing:
            self.explorer._apply_thresholds(self)


class LevelExplorer:
    """An interactive editor of one level's conditions and candidates (see the module).

    level is the Overlaps edited; lengths is (min_len, max_len) of the segments (or
    fragments) its overlaps join, which candidate generation and the gaps need; name is
    "segment" or "fragment", for the texts. show, a and b are those of Fragmentor.explore,
    as is row_height. Works on a bare Overlaps as well.

    It displays as the notebook widget it is built around (widget). refresh() re-reads
    everything from the live objects, after changes made outside the explorer.
    """

    def __init__(self, level: Overlaps, lengths: tuple[int, int], name: str = "segment",
                 show: int | str | list[int | str] | None = None, a: int = 0, b: int | None = None,
                 row_height: int = 280):
        self.level = level
        self.min_len, self.max_len = lengths
        self.name = name
        if row_height < MIN_ROW_HEIGHT:
            raise ValueError(f"row_height must be at least {MIN_ROW_HEIGHT} px, got {row_height}")
        self.row_height = int(row_height)
        n_starts = level.grid.shape[0]
        if b is None:
            b = n_starts
        if not 0 <= a < b <= n_starts:
            raise ValueError(f"View range [{a}, {b}) must be non-empty and within [0, {n_starts})")
        self._x_range = [a - 0.5, b - 0.5]
        self._shown = set(self._pick(show))

        self.fig: go.FigureWidget | None = None
        self._blocks: list[int] = []  # the condition index of each block, top to bottom
        self._geo: _Geometry | None = None
        self._minima: dict[int, np.ndarray] = {}
        self._panels: dict[int, _ConditionPanel] = {}
        self._slots: dict[int, dict[int, int | None]] = {}  # condition -> relaxed id -> colour slot
        self._lane_runs: dict[str, tuple] = {}
        self._gaps: list[tuple[int, int]] = []
        self._stale = level.candidates_stale
        self._syncing = False
        self._bulk = False

        self.status = widgets.HTML()
        self.message = widgets.HTML()

        self.checks: list[widgets.Checkbox] = []
        self.checklist = widgets.VBox(layout=widgets.Layout(
            max_height="136px", overflow="hidden auto", border=f"1px solid {GRID}", padding="2px 6px"))
        self.show_all, self.show_none = _button("all", "plot every condition"), _button("none", "plot none")
        self.show_all.on_click(lambda _: self._show_all(True))
        self.show_none.on_click(lambda _: self._show_all(False))
        conditions_box = widgets.VBox(
            [widgets.HTML("<b>Conditions</b>"), self.checklist, widgets.HBox([self.show_all, self.show_none])],
            layout=widgets.Layout(width="250px", margin="0 16px 0 0"))

        self.min_step = _text("min_step", width=130, label=60)
        self.min_step.observe(self._on_min_step, names="value")
        self.regenerate = _button("Regenerate candidates", "generate the candidates afresh from the conditions")
        self.regenerate.on_click(self._on_regenerate)
        refresh = _button("Refresh", "re-read everything, after changes made outside the explorer")
        refresh.on_click(lambda _: self.refresh())
        self.new_start = _text("start", placeholder="", width=120, label=34)
        self.new_end = _text("end", placeholder="", width=110, label=26)
        self.new_penalty = _text("penalty", placeholder="0", width=120, label=50)
        self.add_candidate = _button("Add candidate", "add this overlap by hand: it survives every regeneration")
        self.add_candidate.on_click(self._on_add_candidate)
        self.manual_list = widgets.VBox()
        self.manual = widgets.Accordion(children=[self.manual_list], selected_index=None)
        candidates_box = widgets.VBox([
            widgets.HTML("<b>Candidates</b>"),
            widgets.HBox([self.min_step, self.regenerate, refresh]),
            widgets.HBox([self.new_start, self.new_end, self.new_penalty, self.add_candidate]),
            self.manual,
        ], layout=widgets.Layout(flex="1 1 auto"))

        self.figure_box = widgets.Box(layout=widgets.Layout(flex="1 1 0%", min_width="0", margin="0"))
        self.panels = widgets.VBox(layout=widgets.Layout(width=f"{PANEL_WIDTH}px", flex="0 0 auto",
                                                         margin="0"))
        self._spacer = widgets.Box(layout=widgets.Layout(height=f"{_Geometry.TOP}px", margin="0"))
        self.widget = widgets.VBox([
            self.status,
            widgets.HBox([conditions_box, candidates_box]),
            self.message,
            widgets.HBox([self.figure_box, self.panels],
                         layout=widgets.Layout(width="100%", align_items="flex-start")),
        ])

        self._sync_checklist()
        self._sync_min_step()
        self._rebuild()

    def _repr_mimebundle_(self, **kwargs):
        return self.widget._repr_mimebundle_(**kwargs)

    def _pick(self, show) -> list[int]:
        """The positions in level.conditions of the conditions show names."""
        conditions = self.level.conditions
        if show is None:
            return list(range(min(DEFAULT_SHOWN, len(conditions))))
        if isinstance(show, (int, str)):
            show = [show]
        picked = []
        for item in show:
            if isinstance(item, str):
                c = self.level.get_condition(item)
                picked.append(next(i for i, other in enumerate(conditions) if other is c))
            elif 0 <= item < len(conditions):
                picked.append(int(item))
            else:
                raise ValueError(f"No condition at position {item}: the level has {len(conditions)}")
        return picked

    def panel(self, index: int) -> _ConditionPanel:
        """The panel of the condition at index in level.conditions, made on first use."""
        if index not in self._panels:
            self._panels[index] = _ConditionPanel(self, index, self.level.conditions[index], self.row_height)
        return self._panels[index]

    def say(self, error: str | None = None, warned=(), info: str | None = None) -> None:
        self.message.value = _notice(error, warned, info)

    def refresh(self) -> None:
        """Re-read everything from the live objects, after changes made outside the explorer:
        new conditions join the checklist, the panels and min_step show the current values,
        and the figure is drawn again."""
        self._sync_checklist()
        self._sync_min_step()
        for panel in self._panels.values():
            panel.sync()
        self._rebuild()
        self.say(info="refreshed")

    # ---- drawing -----------------------------------------------------------------

    def _rebuild(self) -> None:
        """Draw the figure afresh for the ticked conditions, keeping the current zoom."""
        if self.fig is not None and self.fig.layout.xaxis.range is not None:
            self._x_range = list(self.fig.layout.xaxis.range)
        self._blocks = sorted(self._shown)
        self._read_candidates()
        if not self._blocks:
            self.fig, self._geo = None, None
            self.figure_box.children = [widgets.HTML("<i>Tick a condition to plot it.</i>")]
            self.panels.children = []
            return

        self._geo = geo = _Geometry(len(self._blocks), self.row_height)
        layout = dict(height=geo.height, template="none", paper_bgcolor=SURFACE, plot_bgcolor=SURFACE,
                      margin=dict(t=geo.TOP, b=geo.BOTTOM, l=MARGIN_LEFT, r=MARGIN_RIGHT, pad=0,
                                  autoexpand=False),
                      font=dict(family=FONT, size=11, color=INK2), hovermode="closest",
                      dragmode="zoom", showlegend=True, hoverlabel=dict(font_size=11))
        for k in range(geo.n):
            layout.update(_block_layout(geo, k, self.name))
        layout["xaxis"]["range"] = self._x_range

        fig = go.FigureWidget(layout=layout)
        traces = []
        for k, i in enumerate(self._blocks):
            traces.append(self._score_trace(k, i))
            traces += self._light_traces(k, i)
            traces += self._lane_traces(k, i)
        with fig.batch_update():
            fig.add_traces(traces)
            fig.layout.shapes = self._shapes()
            fig.layout.annotations = self._annotations()
        self.fig = fig
        self.figure_box.children = [fig]
        self.panels.children = [self._spacer] + [self.panel(i).widget for i in self._blocks]

    def _score_trace(self, k: int, i: int):
        if i not in self._minima:  # the maps never change, so neither do the minima
            self._minima[i] = self.level.conditions[i].row_minima()
        minima = self._minima[i]
        sx, sy, _, _ = _axes(k)
        return go.Scattergl(
            x=np.arange(len(minima), dtype=np.int32), y=minima, mode="markers", name="min score",
            marker=dict(size=4, color=SCORE), xaxis=sx, yaxis=sy, legend=_legend(k),
            hovertemplate="start %{x:,}<br>min score %{y:.4g}<extra></extra>",
            meta=dict(block=i, kind="score"),
        )

    def _light_traces(self, k: int, i: int) -> list:
        """The threshold lines of condition i and the bounds of its relaxed versions."""
        c = self.level.conditions[i]
        sx, sy, _, _ = _axes(k)
        common = dict(xaxis=sx, yaxis=sy, legend=_legend(k), hoverinfo="skip",
                      meta=dict(block=i, kind="light"))
        traces = []
        x = [-0.5, self.level.grid.shape[0] - 0.5]
        for name, value, color, dash in (("min_val", c.min_val, INK2, "dash"),
                                         ("max_val", c.max_val, INK2, "dashdot"),
                                         ("opt_val", c.opt_val, OPT, "dot")):
            if value is not None:
                traces.append(go.Scatter(x=x, y=[float(value)] * 2, mode="lines",
                                         name=f"{name} = {format_value(value)}",
                                         line=dict(color=color, width=2, dash=dash), **common))
        colors = self._relaxed_colors(i, c)
        for r in sorted(c.relaxed, key=lambda r: r.id):
            rows = relaxed_rows(c, r)
            bounds = [float(v) for v in (r.min_val, r.max_val) if v is not None]
            if rows is None or not bounds:
                continue
            x0, x1 = rows[0] - 0.5, rows[1] + 0.5
            xs, ys, labels = [], [], []
            for v in bounds:  # each bound a segment over the starts it applies to, labelled at its end
                xs += [x0, x1, None]
                ys += [v, v, None]
                labels += ["", f"#{r.id}", ""]
            traces.append(go.Scatter(
                x=xs, y=ys, mode="lines+text", text=labels, textposition="top left",
                textfont=dict(size=10, color=INK2), line=dict(color=colors[r.id], width=2),
                name=f"relaxed {describe_relaxed(r)}", **common))
        return traces

    def _read_candidates(self) -> None:
        """Split the level's candidates into the lanes, find the gaps, note staleness."""
        level = self.level
        width = len(level.seq) + 1  # (start, end) -> one key
        m = level.manual
        m_starts, m_ends = m.starts[1:-1], m.ends[1:-1]
        runs = {"manual": coverage(m_starts, m_ends)}
        c = level.candidates
        if c is None:
            runs["pass"] = runs["relaxed"] = coverage([], [])
            self._gaps = []
        else:
            starts, ends, penalty = c.starts[1:-1], c.ends[1:-1], c.extra_penalty[1:-1]
            generated = ~np.isin(starts * width + ends, m_starts * width + m_ends)
            runs["pass"] = coverage(starts[generated & (penalty == 0)], ends[generated & (penalty == 0)])
            runs["relaxed"] = coverage(starts[generated & (penalty != 0)], ends[generated & (penalty != 0)])
            self._gaps = level.get_gapped_regions(self.max_len + 1)
        self._lane_runs = runs
        self._stale = level.candidates_stale
        self._update_status()
        self._sync_manual()

    def _lane_traces(self, k: int, i: int) -> list:
        """The candidate lanes of block k, each run of overlapping candidates one thick
        stroke, and the hover targets of the gaps."""
        _, _, lx, ly = _axes(k)
        common = dict(xaxis=lx, yaxis=ly, showlegend=False, hovertemplate="%{customdata}<extra></extra>")
        faded = self._stale and self.level.candidates is not None
        traces = []
        for y, lane in zip((2, 1, 0), LANES):
            starts, ends, counts = self._lane_runs[lane]
            if not starts.size:
                continue
            nan = np.full(starts.size, np.nan)
            xs = np.column_stack([starts, (starts + ends) / 2, ends, nan]).ravel()
            labels = [f"{lane}: [{s:,}, {e:,}) · {n:,} candidate{'' if n == 1 else 's'}"
                      for s, e, n in zip(starts.tolist(), ends.tolist(), counts.tolist())]
            hover = np.array([[t, t, t, None] for t in labels], dtype=object).ravel()
            traces.append(go.Scatter(x=xs, y=np.where(np.isnan(xs), np.nan, y), mode="lines",
                                     line=dict(color=MUTED if faded else INK2, width=8),
                                     opacity=0.55 if faded else 1.0, customdata=hover,
                                     meta=dict(block=i, kind="lanes"), **common))
        if self._gaps:
            traces.append(go.Scatter(
                x=[(a + b) / 2 for a, b in self._gaps], y=[1] * len(self._gaps), mode="markers",
                marker=dict(size=16, opacity=0),
                customdata=[f"gap [{a:,}, {b:,}): crossing it takes a {self.name} of {b - a:,} bp, "
                            f"longer than {self.max_len:,}" for a, b in self._gaps],
                meta=dict(block=i, kind="lanes"), **common))
        return traces

    def _shapes(self) -> list[dict]:
        """The area washes of the blocks and the gap washes in their lanes."""
        shapes = []
        for k, i in enumerate(self._blocks):
            sx, sy, lx, ly = _axes(k)
            rows = self.level.conditions[i].area_rows()
            if rows is not None:
                shapes.append(dict(type="rect", xref=sx, yref=f"{sy} domain", x0=rows[0] - 0.5,
                                   x1=rows[1] + 0.5, y0=0, y1=1, fillcolor=AREA, opacity=0.10,
                                   line_width=0, layer="below"))
            for a, b in self._gaps:
                shapes.append(dict(type="rect", xref=lx, yref=f"{ly} domain", x0=a, x1=b, y0=0, y1=1,
                                   fillcolor=CRITICAL, opacity=0.22, line_width=0, layer="below"))
        return shapes

    def _annotations(self) -> list[dict]:
        """The block titles, and a note on the lanes when the candidates are missing or stale."""
        if self.level.candidates is None:
            note = "no candidates yet: press Regenerate"
        elif self._stale:
            note = "stale: the conditions changed since, press Regenerate"
        else:
            note = None
        out = []
        for k, i in enumerate(self._blocks):
            c = self.level.conditions[i]
            area = "" if c.area_rows() is None else f"  <span style='color:{MUTED}'>area [{format_value(c.start) or '…'}, {format_value(c.end) or '…'})</span>"
            out.append(dict(text=f"<b>[{i}] {html.escape(c.label)}</b>{area}", x=0, xref="paper",
                            xanchor="left", y=self._geo.top(k), yref="paper", yanchor="top",
                            showarrow=False, font=dict(size=12, color=INK)))
            if note:
                _, _, lx, ly = _axes(k)
                out.append(dict(text=note, x=1, xref=f"{lx} domain", xanchor="right", y=1,
                                yref=f"{ly} domain", yanchor="bottom", showarrow=False,
                                font=dict(size=10, color=INK2), bgcolor="rgba(252,252,251,0.85)"))
        return out

    def _relayout(self) -> None:
        self.fig.layout.shapes = self._shapes()
        self.fig.layout.annotations = self._annotations()

    def _redraw_block(self, i: int) -> None:
        """Draw condition i's thresholds, relaxed bounds and area again."""
        if self.fig is None or i not in self._blocks:
            return
        new = self._light_traces(self._blocks.index(i), i)
        with self.fig.batch_update():
            self.fig.data = [t for t in self.fig.data if (t.meta["block"], t.meta["kind"]) != (i, "light")]
            if new:
                self.fig.add_traces(new)
            self._relayout()

    def _redraw_candidates(self) -> None:
        """Read the candidates again and draw every block's lanes again."""
        self._read_candidates()
        if self.fig is None:
            return
        new = [t for k, i in enumerate(self._blocks) for t in self._lane_traces(k, i)]
        with self.fig.batch_update():
            self.fig.data = [t for t in self.fig.data if t.meta["kind"] != "lanes"]
            if new:
                self.fig.add_traces(new)
            self._relayout()

    def _changed(self, panel: _ConditionPanel) -> None:
        """After an edit of panel's condition: its block, and the lanes once they go stale."""
        self._redraw_block(panel.index)
        if self.level.candidates_stale != self._stale:
            self._redraw_candidates()
        else:
            self._update_status()

    def _relaxed_colors(self, i: int, c: Condition) -> dict[int, str]:
        """The colour of each relaxed version of condition i. A version takes the first free
        slot when the explorer first meets it and keeps it as long as it exists."""
        slots = self._slots.setdefault(i, {})
        current = {r.id for r in c.relaxed}
        for rid in list(slots):
            if rid not in current:
                del slots[rid]
        for rid in sorted(current - slots.keys()):
            used = set(slots.values())
            slots[rid] = next((s for s in range(len(RELAXED)) if s not in used), None)
        return {rid: MUTED if s is None else RELAXED[s] for rid, s in slots.items()}

    # ---- texts ---------------------------------------------------------------------

    def _update_status(self) -> None:
        level = self.level
        head = (f"<b>{self.name.capitalize()} level</b> · {len(level.seq):,} bp · overlaps "
                f"{level.min_len}–{level.max_len} bp · {self.name}s {self.min_len:,}–{self.max_len:,} bp "
                f"· {len(level.conditions)} condition{'' if len(level.conditions) == 1 else 's'}")
        c = level.candidates
        if c is None:
            line = "no candidates generated yet"
        else:
            penalised = int((c.extra_penalty[1:-1] != 0).sum())
            line = (f"{c.n_overlaps:,} candidates ({penalised:,} with a penalty, "
                    f"{level.manual.n_overlaps:,} added by hand)")
            if self._gaps:
                shown = ", ".join(f"[{a:,}, {b:,})" for a, b in self._gaps[:5])
                more = " …" if len(self._gaps) > 5 else ""
                line += (f" · <span style='color:{CRITICAL}'>✖ {len(self._gaps)} "
                         f"gap{'' if len(self._gaps) == 1 else 's'} no {self.name} can cross</span>: "
                         f"{shown}{more}")
            else:
                line += f" · no gap is too long for a {self.name}"
            if level.candidates_stale:
                line += " · <b>stale</b>: the conditions changed since, press Regenerate"
        self.status.value = f"{head}<br>{line}"

    def _sync_checklist(self) -> None:
        conditions = self.level.conditions
        for i in range(len(self.checks), len(conditions)):  # conditions added since
            box = widgets.Checkbox(value=i in self._shown, indent=False,
                                   layout=widgets.Layout(width="auto", margin="0"))
            box.observe(lambda change, i=i: self._on_check(i, change["new"]), names="value")
            self.checks.append(box)
        for i, box in enumerate(self.checks):
            c = conditions[i]
            box.description = f"[{i}] {c.label}" + ("" if c.is_constraint else " (objective)")
        self.checklist.children = self.checks or [widgets.HTML("<i>no conditions on this level</i>")]

    def _sync_min_step(self) -> None:
        self._syncing = True
        try:
            self.min_step.value = format_value(self.level.min_step)
        finally:
            self._syncing = False

    def _sync_manual(self) -> None:
        m = self.level.manual
        pairs = list(zip(m.starts[1:-1].tolist(), m.ends[1:-1].tolist(), m.extra_penalty[1:-1].tolist()))
        rows = []
        for s, e, p in pairs[:MAX_LISTED]:
            remove = _button("✕", f"remove ({s}, {e})", width="32px")
            remove.on_click(lambda _, s=s, e=e: self._remove_manual(s, e))
            rows.append(widgets.HBox([remove, widgets.HTML(f"({s:,}, {e:,}) · {e - s} bp · +{p:g}")]))
        if len(pairs) > MAX_LISTED:
            rows.append(widgets.HTML(f"… and {len(pairs) - MAX_LISTED:,} more "
                                     f"(Fragmentor.remove_candidates takes any of them)"))
        self.manual_list.children = rows or [widgets.HTML("<i>none</i>")]
        self.manual.set_title(0, f"Added by hand ({len(pairs)})")

    # ---- actions -------------------------------------------------------------------

    def _on_check(self, i: int, on: bool) -> None:
        (self._shown.add if on else self._shown.discard)(i)
        if not self._bulk:
            self._rebuild()

    def _show_all(self, on: bool) -> None:
        self._bulk = True
        try:
            for box in self.checks:
                box.value = on
        finally:
            self._bulk = False
        self._rebuild()

    def _apply_thresholds(self, panel: _ConditionPanel) -> None:
        """Send the whole panel through one Condition.update, so that a pair of changes
        that only works together (min_val above the old max_val) goes through once both
        are typed. A refused set changes nothing and leaves what was typed."""
        try:
            values = {name: parse_value(f.value, int if name in ("start", "end") else float)
                      for name, f in panel.fields.items()}
        except ValueError as e:
            panel.say(error=str(e))
            return
        error, warned = _attempt(lambda: panel.condition.update(**values))
        if error:
            panel.say(error=error, warned=warned)
            return
        panel.sync()
        panel.say(warned=warned)
        self._changed(panel)

    def _add_relaxed(self, panel: _ConditionPanel) -> None:
        new = panel.new
        try:
            kwargs = dict(min_val=parse_value(new["min_val"].value), max_val=parse_value(new["max_val"].value),
                          extra_penalty=parse_value(new["penalty"].value) or 0.0,
                          start=parse_value(new["start"].value, int), end=parse_value(new["end"].value, int))
        except ValueError as e:
            panel.say(error=str(e))
            return
        added = []
        error, warned = _attempt(lambda: added.append(panel.condition.add_relaxed(**kwargs)))
        if error:
            panel.say(error=error, warned=warned)
            return
        for f in new.values():
            f.value = ""
        panel.sync_relaxed()
        panel.say(warned=warned, info=f"added relaxed #{added[0]}")
        self._changed(panel)

    def _remove_relaxed(self, panel: _ConditionPanel, rid: int) -> None:
        error, warned = _attempt(lambda: panel.condition.remove_relaxed(rid))
        if error:
            panel.say(error=error, warned=warned)
            return
        panel.sync_relaxed()
        panel.say(warned=warned, info=f"removed relaxed #{rid}")
        self._changed(panel)

    def _on_min_step(self, change) -> None:
        if self._syncing:
            return
        try:
            value = parse_value(self.min_step.value, int)
        except ValueError as e:
            self.say(error=str(e))
            return
        error, warned = _attempt(lambda: setattr(self.level, "min_step", value))
        if error:
            self.say(error=error, warned=warned)
            return
        self._sync_min_step()  # a step of 1 is stored as None
        step = self.level.min_step
        self.say(warned=warned, info=f"min_step is now {'none (no thinning)' if step is None else step}"
                                     + ("; press Regenerate to apply it" if self.level.candidates_stale else ""))
        self._redraw_candidates()

    def _on_regenerate(self, _) -> None:
        self.regenerate.disabled = True
        self.say(info="regenerating the candidates…")
        began = time.perf_counter()
        try:
            error, warned = _attempt(lambda: self.level.generate_candidates(self.min_len))
        finally:
            self.regenerate.disabled = False
        if error:
            self.say(error=error, warned=warned)
            return
        self.say(warned=warned, info=f"generated {self.level.candidates.n_overlaps:,} candidates "
                                     f"in {time.perf_counter() - began:.2f} s")
        self._redraw_candidates()

    def _on_add_candidate(self, _) -> None:
        try:
            start, end = parse_value(self.new_start.value, int), parse_value(self.new_end.value, int)
            penalty = parse_value(self.new_penalty.value)
        except ValueError as e:
            self.say(error=str(e))
            return
        if start is None or end is None:
            self.say(error="a candidate needs a start and an end")
            return
        penalty = 0.0 if penalty is None else penalty
        if not penalty >= 0:  # also refuses NaN, as for relaxed versions
            self.say(error=f"penalty {penalty:g} must be non-negative")
            return
        error, warned = _attempt(lambda: self.level.add_candidates(start, end, penalty))
        if error:
            self.say(error=error, warned=warned)
            return
        lo, hi = self.level.min_len, self.level.max_len
        if not lo <= end - start <= hi:
            warned.append(f"({start:,}, {end:,}) is {end - start} bp long, outside this level's "
                          f"overlap lengths {lo}–{hi}; added anyway")
        self.new_start.value = self.new_end.value = ""
        self.say(warned=warned, info=f"added ({start:,}, {end:,}) by hand, penalty {penalty:g}")
        self._redraw_candidates()

    def _remove_manual(self, start: int, end: int) -> None:
        error, warned = _attempt(lambda: self.level.remove_candidates(start, end))
        if error:
            self.say(error=error, warned=warned)
            return
        self.say(warned=warned, info=f"removed ({start:,}, {end:,})")
        self._redraw_candidates()
