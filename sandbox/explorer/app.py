"""Interactive explorer for choppy assembly-split results.

Streamlit app to explore the output of the local-neighbourhood splitting
algorithm (see ``sandbox/nopool_localneigh.ipynb``). It loads two pickles:

* ``full_all_overlaps.pkl``  -- {seq_id: [overlap, ...]} position-sorted nodes.
      ``seg`` overlaps are segment junctions (cross-sequence homology,
      range == (0, seq_len)); ``fr`` overlaps are fragment overlaps with
      primers and a *local* homology-free ``range``. Two sentinels bookend
      the list: (0, 0) and (seq_len, seq_len).
* ``full_overlap_states.pkl`` -- {seq_id: [state_dict, ...]} indexed the same
      way as all_overlaps. ``state_dict`` maps
      ``state = (seg_start, range_end)`` -> ``(weight, pred_index, pred_state)``.
      ``seg_start`` is where the current segment began; ``range_end`` is the
      running-min right edge of the local neighbourhood. The optimal path
      ending at a given (index, state) is recovered by walking predecessors
      back to pred_index == -1.

Run with::

    .venv/bin/streamlit run sandbox/explorer/app.py

(from the repo root, or set CHOPPY_DATA_DIR / CHOPPY_FASTA).
"""
from __future__ import annotations

import os
import pickle
from pathlib import Path

import plotly.graph_objects as go
import streamlit as st

# --------------------------------------------------------------------------
# Config / paths
# --------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("CHOPPY_DATA_DIR", REPO_ROOT / "sandbox" / "data"))
FASTA = Path(os.environ.get("CHOPPY_FASTA", DATA_DIR / "20240414-forSveta.fa"))
ALL_OVERLAPS_PKL = DATA_DIR / "full_all_overlaps.pkl"
OVERLAP_STATES_PKL = DATA_DIR / "full_overlap_states.pkl"

# Below this window width (bp) we render individual bases and primer bars.
BASE_THRESHOLD = 250

# Lane y-centres for the genomic track.
LANES = {
    "bases":     5.0,   # sequence ruler / A C G T letters when zoomed
    "seg":       4.0,   # candidate seg overlaps in range
    "fr":        3.0,   # candidate fr overlaps in range
    "primer":    2.4,   # primers (selected overlap, or all when zoomed)
    "path_seg":  1.5,   # segments of the reconstructed path
    "path_frag": 0.9,   # fragments (edges) of the reconstructed path
    "path_node": 0.3,   # overlap nodes used by the path
}

COL_SEG = "#3273dc"       # blue  - seg overlaps
COL_FR = "#ff9f1c"        # orange - fr overlaps
COL_SELECTED = "#e63946"  # red   - selected overlap
COL_SEGMENT = "#8ecae6"   # segment block
COL_FRAGMENT = "#95d5b2"  # fragment block
COL_RANGE_SEL = "rgba(230, 57, 70, 0.08)"     # selected view range
COL_RANGE_NBR = "rgba(46, 196, 182, 0.18)"    # local-neighbourhood range


# --------------------------------------------------------------------------
# Data loading (cached once per process)
# --------------------------------------------------------------------------
@st.cache_resource(show_spinner="Loading overlaps + states (this is large, ~1 min)…")
def load_data():
    with open(ALL_OVERLAPS_PKL, "rb") as fh:
        all_overlaps = pickle.load(fh)
    with open(OVERLAP_STATES_PKL, "rb") as fh:
        overlap_states = pickle.load(fh)

    sequences = {}
    if FASTA.exists():
        from Bio import SeqIO
        for rec in SeqIO.parse(str(FASTA), "fasta"):
            sequences[rec.id] = str(rec.seq).upper()
    return all_overlaps, overlap_states, sequences


def seq_length(all_overlaps, sequences, seq_id):
    """Sequence length from the FASTA if present, else the end sentinel."""
    if seq_id in sequences:
        return len(sequences[seq_id])
    return max(ov["pos"][1] for ov in all_overlaps[seq_id])


# --------------------------------------------------------------------------
# Path reconstruction
# --------------------------------------------------------------------------
def reconstruct_path(states, index, state):
    """Walk predecessors from (index, state) back to the start sentinel.

    Returns a list of (overlap_index, state) ordered start -> end.
    """
    path = []
    cur_i, cur_s = index, state
    seen = set()
    while cur_i != -1:
        if (cur_i, cur_s) in seen:      # guard against malformed cycles
            break
        seen.add((cur_i, cur_s))
        path.append((cur_i, cur_s))
        _, pred_i, pred_s = states[cur_i][cur_s]
        cur_i, cur_s = pred_i, pred_s
    return path[::-1]


def segment_blocks(path, overlaps):
    """Group a path into segments.

    A node is a segment boundary when its seg_start (state[0]) equals the
    overlap's own start position -- mirroring the notebook's annotator.
    Returns list of dicts: {start, end, nodes:[path_idx...]}.
    """
    boundaries = [
        pi for pi, (ov_i, state) in enumerate(path)
        if state[0] == overlaps[ov_i]["pos"][0]
    ]
    if path and boundaries and boundaries[-1] != len(path) - 1:
        boundaries.append(len(path) - 1)  # partial path ending mid-segment

    segs = []
    for si in range(len(boundaries) - 1):
        pa, pb = boundaries[si], boundaries[si + 1]
        left = overlaps[path[pa][0]]
        right = overlaps[path[pb][0]]
        segs.append({
            "start": max(left["pos"][0], 0),
            "end": right["pos"][1],
            "nodes": list(range(pa, pb + 1)),
        })
    return segs


# --------------------------------------------------------------------------
# Figure building
# --------------------------------------------------------------------------
def _hline_traces(spans, yc, color, name, hovertexts, width=6, height=0.32):
    """Two traces: line segments for the bars + hover markers at midpoints."""
    xs, ys = [], []
    mx, my, mtext = [], [], []
    for (x0, x1), ht in zip(spans, hovertexts):
        xs += [x0, x1, None]
        ys += [yc, yc, None]
        mx.append((x0 + x1) / 2)
        my.append(yc)
        mtext.append(ht)
    line = go.Scattergl(
        x=xs, y=ys, mode="lines", line=dict(color=color, width=width),
        name=name, hoverinfo="skip", showlegend=True,
    )
    markers = go.Scattergl(
        x=mx, y=my, mode="markers", marker=dict(color=color, size=8, opacity=0.35),
        name=name, hoverinfo="text", hovertext=mtext, showlegend=False,
    )
    return [line, markers]


def _rect(fig, x0, x1, yc, color, height=0.34, line_color=None, opacity=0.9):
    fig.add_shape(
        type="rect", x0=x0, x1=x1, y0=yc - height, y1=yc + height,
        fillcolor=color, opacity=opacity,
        line=dict(color=line_color or color, width=1),
        layer="below",
    )


def fr_hover(ov, i):
    return (f"[{i}] fr {ov['pos']}<br>range {ov['range']}"
            f"<br>fwd {ov['pr_left_seq']} (Tm {ov['tm_left']:.1f})"
            f"<br>rev {ov['pr_right_seq']} (Tm {ov['tm_right']:.1f})")


def seg_hover(ov, i):
    return f"[{i}] seg {ov['pos']}<br>range {ov['range']}"


def build_figure(seq_id, overlaps, seq_str, view, in_range, selected_i,
                 path, segs, sel_state, show_candidates):
    v0, v1 = view
    fig = go.Figure()

    # --- selected view range shading ---
    fig.add_vrect(x0=v0, x1=v1, fillcolor=COL_RANGE_SEL, line_width=0, layer="below")

    # --- local-neighbourhood range of the selected state ---
    if sel_state is not None:
        ns, ne = sel_state[0], sel_state[1]
        fig.add_vrect(
            x0=ns, x1=ne, fillcolor=COL_RANGE_NBR, line_width=0, layer="below",
            annotation_text=f"local neighbourhood [{ns}, {ne}]",
            annotation_position="top left",
        )

    # --- candidate overlaps within range ---
    if show_candidates:
        seg_spans, seg_hovers, fr_spans, fr_hovers = [], [], [], []
        for i, ov in in_range:
            if ov["type"] == "seg":
                seg_spans.append(ov["pos"]); seg_hovers.append(seg_hover(ov, i))
            else:
                fr_spans.append(ov["pos"]); fr_hovers.append(fr_hover(ov, i))
        if seg_spans:
            for t in _hline_traces(seg_spans, LANES["seg"], COL_SEG,
                                   "seg overlaps", seg_hovers):
                fig.add_trace(t)
        if fr_spans:
            for t in _hline_traces(fr_spans, LANES["fr"], COL_FR,
                                   "fr overlaps", fr_hovers):
                fig.add_trace(t)

    # --- reconstructed path: segments ---
    for k, s in enumerate(segs):
        _rect(fig, s["start"], s["end"], LANES["path_seg"],
              COL_SEGMENT, height=0.30, line_color="#219ebc", opacity=0.55)
        fig.add_annotation(x=(s["start"] + s["end"]) / 2, y=LANES["path_seg"],
                           text=f"seg{k + 1}", showarrow=False,
                           font=dict(size=10, color="#023047"))

    # --- reconstructed path: fragments (edges) ---
    for j in range(1, len(path)):
        l_ov = overlaps[path[j - 1][0]]
        r_ov = overlaps[path[j][0]]
        fs, fe = max(l_ov["pos"][0], 0), r_ov["pos"][1]
        if fs >= fe:
            continue
        off = 0.16 if j % 2 else -0.16   # stagger so overlaps between frags show
        _rect(fig, fs, fe, LANES["path_frag"] + off,
              COL_FRAGMENT, height=0.14, line_color="#40916c", opacity=0.8)

    # --- reconstructed path: overlap nodes ---
    nx, ny, ntext, ncol, nsym = [], [], [], [], []
    for ov_i, state in path:
        ov = overlaps[ov_i]
        nx.append((ov["pos"][0] + ov["pos"][1]) / 2)
        ny.append(LANES["path_node"])
        is_boundary = state[0] == ov["pos"][0]
        ntext.append(f"[{ov_i}] {ov['type']} {ov['pos']}<br>state {state}"
                     f"<br>{'segment start' if is_boundary else 'interior'}")
        ncol.append(COL_SEG if ov["type"] == "seg" else COL_FR)
        nsym.append("star" if is_boundary else "circle")
    if nx:
        fig.add_trace(go.Scattergl(
            x=nx, y=ny, mode="markers",
            marker=dict(color=ncol, size=11, symbol=nsym,
                        line=dict(color="#333", width=1)),
            name="path nodes", hoverinfo="text", hovertext=ntext, showlegend=True,
        ))

    # --- selected overlap: highlight + primers ---
    if selected_i is not None:
        ov = overlaps[selected_i]
        yc = LANES["seg"] if ov["type"] == "seg" else LANES["fr"]
        _rect(fig, ov["pos"][0], ov["pos"][1], yc, COL_SELECTED,
              height=0.40, line_color=COL_SELECTED, opacity=0.35)
        if ov["type"] == "fr":
            p0, p1 = ov["pos"]
            fwd = (p0, p0 + ov["pr_left_len"])
            rev = (p1 - ov["pr_right_len"], p1)
            for span, col, lbl in [(fwd, "#2a9d8f", "fwd"), (rev, "#9d4edd", "rev")]:
                _rect(fig, span[0], span[1], LANES["primer"], col,
                      height=0.22, opacity=0.9)
            fig.add_annotation(x=(p0 + p1) / 2, y=LANES["primer"],
                               text="primers", showarrow=False,
                               font=dict(size=9, color="#555"))

    # --- bases (only when zoomed in) ---
    width = v1 - v0
    if seq_str and width <= BASE_THRESHOLD:
        lo, hi = max(0, int(v0)), min(len(seq_str), int(v1) + 1)
        xs = list(range(lo, hi))
        letters = list(seq_str[lo:hi])
        fig.add_trace(go.Scattergl(
            x=[x + 0.5 for x in xs], y=[LANES["bases"]] * len(xs),
            mode="text", text=letters, textfont=dict(size=11, family="monospace"),
            hoverinfo="skip", showlegend=False,
        ))
    else:
        fig.add_annotation(x=(v0 + v1) / 2, y=LANES["bases"],
                           text=f"◀ zoom in below {BASE_THRESHOLD} bp to see bases ▶",
                           showarrow=False, font=dict(size=11, color="#999"))

    # --- layout ---
    fig.update_xaxes(range=[v0, v1], title="position (bp)")
    fig.update_yaxes(
        range=[-0.3, 5.6],
        tickmode="array",
        tickvals=[LANES[k] for k in
                  ["bases", "seg", "fr", "primer", "path_seg", "path_frag", "path_node"]],
        ticktext=["sequence", "seg cand.", "fr cand.", "primers",
                  "path segments", "path fragments", "path nodes"],
    )
    fig.update_layout(
        height=560, margin=dict(l=110, r=20, t=30, b=40),
        legend=dict(orientation="h", yanchor="bottom", y=1.01, x=0),
        hovermode="closest", dragmode="pan",
    )
    return fig


# --------------------------------------------------------------------------
# UI
# --------------------------------------------------------------------------
def main():
    st.set_page_config(page_title="Choppy split explorer", layout="wide")
    st.title("Choppy — assembly split explorer")

    for pkl in (ALL_OVERLAPS_PKL, OVERLAP_STATES_PKL):
        if not pkl.exists():
            st.error(f"Missing data file: {pkl}\n\n"
                     "Set CHOPPY_DATA_DIR to the folder holding the pkl files.")
            st.stop()

    all_overlaps, overlap_states, sequences = load_data()

    # ---- sidebar controls ----
    sb = st.sidebar
    sb.header("Selection")
    seq_id = sb.selectbox("Sequence", sorted(all_overlaps.keys()))
    overlaps = all_overlaps[seq_id]
    states = overlap_states[seq_id]
    seq_str = sequences.get(seq_id, "")
    slen = seq_length(all_overlaps, sequences, seq_id)

    sb.caption(f"{len(overlaps):,} overlap nodes · length {slen:,} bp"
               + ("" if seq_str else " · (FASTA not found — no bases)"))

    default_end = min(2000, slen)
    view = sb.slider("Range / view window (bp)", 0, slen, (0, default_end), step=1)
    v0, v1 = view
    sb.caption(f"window = {v1 - v0:,} bp"
               + ("  ·  bases shown" if (v1 - v0) <= BASE_THRESHOLD and seq_str
                  else f"  ·  zoom < {BASE_THRESHOLD} bp for bases"))

    show_candidates = sb.checkbox("Show candidate overlaps in range", value=True)
    type_filter = sb.multiselect("Overlap types", ["seg", "fr"], default=["seg", "fr"])

    # overlaps whose position falls within the selected range
    in_range = [
        (i, ov) for i, ov in enumerate(overlaps)
        if ov["pos"][0] >= v0 and ov["pos"][1] <= v1 and ov["type"] in type_filter
    ]
    sb.caption(f"{len(in_range):,} overlap(s) in range")

    # ---- pick an overlap in range ----
    sb.header("Overlap → states → path")
    if in_range:
        labels = [
            f"[{i}] {ov['type']} ({ov['pos'][0]},{ov['pos'][1]})"
            + (f"  ·states {len(states[i])}" if states[i] else "  ·no states")
            for i, ov in in_range
        ]
        # default to the first overlap that actually has states
        default_idx = next((k for k, (i, _) in enumerate(in_range) if states[i]), 0)
        pick = sb.selectbox("Overlap in range", range(len(in_range)),
                            index=default_idx, format_func=lambda k: labels[k])
        selected_i, selected_ov = in_range[pick]
    else:
        selected_i, selected_ov = None, None
        sb.info("No overlaps in the selected range.")

    # ---- pick a state of that overlap ----
    sel_state = None
    path, segs = [], []
    if selected_i is not None:
        st_dict = states[selected_i]
        if st_dict:
            # sort states by weight ascending (best first)
            st_items = sorted(st_dict.items(), key=lambda kv: kv[1][0])
            st_labels = [
                f"seg_start={s[0]}, range_end={s[1]}  ·  w={v[0]:.1f}"
                for s, v in st_items
            ]
            sk = sb.selectbox("State (of selected overlap)", range(len(st_items)),
                             format_func=lambda k: st_labels[k])
            sel_state = st_items[sk][0]
            path = reconstruct_path(states, selected_i, sel_state)
            segs = segment_blocks(path, overlaps)
        else:
            sb.warning("This overlap is unreachable (no states) — no path to draw.")

    # ---- main figure ----
    fig = build_figure(seq_id, overlaps, seq_str, view, in_range, selected_i,
                       path, segs, sel_state, show_candidates)
    st.plotly_chart(fig, use_container_width=True,
                    config={"scrollZoom": True, "displaylogo": False})

    st.caption("Drag to pan · scroll to zoom · shaded red = selected range · "
               "shaded teal = local neighbourhood of the chosen state · "
               "★ path node = segment start.")

    # ---- details ----
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("Selected overlap")
        if selected_ov is not None:
            st.write(f"index **{selected_i}** · type **{selected_ov['type']}** · "
                     f"pos **{selected_ov['pos']}** · range **{selected_ov['range']}**")
            if selected_ov["type"] == "fr":
                st.write(f"fwd primer `{selected_ov['pr_left_seq']}` "
                         f"(len {selected_ov['pr_left_len']}, Tm {selected_ov['tm_left']:.1f})")
                st.write(f"rev primer `{selected_ov['pr_right_seq']}` "
                         f"(len {selected_ov['pr_right_len']}, Tm {selected_ov['tm_right']:.1f})")
                st.write(f"boundary_motif: {selected_ov.get('boundary_motif')}")
        else:
            st.write("—")

    with c2:
        st.subheader("Chosen state")
        if sel_state is not None:
            w, pi, ps = states[selected_i][sel_state]
            st.write(f"state **(seg_start={sel_state[0]}, range_end={sel_state[1]})**")
            st.write(f"weight **{w:.1f}** · predecessor index **{pi}** · "
                     f"predecessor state **{ps}**")
            st.write(f"local neighbourhood: **[{sel_state[0]}, {sel_state[1]}]** "
                     f"({sel_state[1] - sel_state[0]:,} bp)")
        else:
            st.write("—")

    # ---- path table ----
    if path:
        st.subheader(f"Reconstructed path — {len(segs)} segment(s), "
                     f"{max(len(path) - 1, 0)} fragment(s)")
        rows = []
        for step, (ov_i, state) in enumerate(path):
            ov = overlaps[ov_i]
            w = states[ov_i][state][0]
            rows.append({
                "step": step,
                "ov_index": ov_i,
                "type": ov["type"],
                "pos": f"{ov['pos'][0]}–{ov['pos'][1]}",
                "seg_start": state[0],
                "range_end": state[1],
                "cum_weight": round(w, 1),
                "segment_start?": "★" if state[0] == ov["pos"][0] else "",
            })
        st.dataframe(rows, use_container_width=True, hide_index=True)


if __name__ == "__main__":
    main()
