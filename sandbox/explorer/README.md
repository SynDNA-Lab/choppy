# Choppy split explorer

Interactive local app to explore the assembly-split results produced by
`sandbox/nopool_localneigh.ipynb`.

## Run

From the repo root:

```bash
.venv/bin/streamlit run sandbox/explorer/app.py
```

Then open the URL it prints (default http://localhost:8501). The first load
unpickles `full_overlap_states.pkl` (~840 MB), which takes ~1 min and stays
cached for the life of the server process.

### Data locations

By default it reads from `sandbox/data/`:

- `full_all_overlaps.pkl` — overlap nodes per sequence
- `full_overlap_states.pkl` — DP states per sequence
- `20240414-forSveta.fa` — sequences (for showing bases; optional)

Override with env vars if needed:

```bash
CHOPPY_DATA_DIR=/path/to/pkls CHOPPY_FASTA=/path/to/seqs.fa \
  .venv/bin/streamlit run sandbox/explorer/app.py
```

## What you can do

- **Pick a sequence** and **set a range / view window** (sidebar slider). Zoom
  below 250 bp and the actual DNA bases appear; the selected `fr` overlap's
  primers are drawn as bars.
- **See all candidate overlaps in range** on the track: `seg` (segment
  junctions, blue) and `fr` (fragment overlaps, orange). Filter by type.
- **Pick an overlap → pick one of its states**, and the app reconstructs and
  draws that state's optimal path from the sequence start: segments, fragments
  (staggered so their shared overlaps are visible) and the overlap nodes
  (★ = segment start).
- The **local neighbourhood** of the chosen state — `[seg_start, range_end]`,
  the window a segment's fragments must lie within — is shaded in teal.
- Details panels + a full path table (per-step cumulative weight) below the
  plot.

Pan by dragging, zoom by scrolling. Re-running a control resets the plot to the
sidebar range.
