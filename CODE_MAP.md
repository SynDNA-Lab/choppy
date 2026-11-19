# Code map: `fragmentor.py` and `constraints.py`

A reference for the classes, methods and helpers in
[`choppy/fragmentor.py`](choppy/fragmentor.py) and [`choppy/constraints.py`](choppy/constraints.py),
private (`_`-prefixed) members included, followed by how the pieces connect.

---

## `choppy/constraints.py`

### Module-level helpers

| Function | What it does |
|---|---|
| `reverse_complement(seq)` | Returns the reverse complement of an ACGT string. |
| `invalid_bases(seq)` | Returns the positions of characters that aren't A, C, G or T. |
| `check_seq(seq, name)` | Raises if the sequence has non-ACGT characters, otherwise returns it unchanged. Used by `Fragmentor.__post_init__`. |
| `valid_kmers(seq, kmer)` | Mask of the k-mers made only of ACGT. Used by `BackgroundKmers.encode_seqs` to drop ambiguous k-mers. |
| `encode_kmers(seq, kmer)` | 2-bit base-4 code for each k-mer (as uint64). Used by `BackgroundKmers.encode_seqs` and `lookup`. |
| `_same(a, b)` | Compares two threshold values, treating `None` correctly. Used by `Condition.update`. |
| `THRESHOLDS`, `_UNSET` | The names of the threshold fields, and the "not given" sentinel for `Condition.update` (where `None` means "remove"). |

### `OverlapGrid` (frozen dataclass)

Converts between overlap coordinates and map cells. Cell `(i, j)` is the overlap
`(i, i + min_len + j)`. Fields: `seq_len`, `min_len`, `max_len`.

- `shape`: `(n_starts, n_lengths)`.
- `check_cell(i, j)`: raises if a cell is out of bounds.
- `cell(start, end)`: converts an overlap to a cell. The end isn't bounded, so cells
  that run past the end of the sequence still exist.
- `overlap(i, j)`: converts a cell back to `(start, end)`.
- `starts` / `lengths` / `ends`: broadcastable arrays covering the whole grid, so a rule
  can be written as one vectorised expression.
- `slices(sub)`: the row and column slices of this grid that hold a narrower sub-grid.
  This is what lets the two levels share maps.
- `edge_mask(min_edge_len)`: drops overlaps too close to either edge for a
  minimum-length segment or fragment to fit.

### `OverlapMap` (dataclass)

One scorer's scores over a grid. Fields: `grid`, `values` (float16 in 0..1),
`multiplier` (turns scores back into the scorer's own units).

- `get_cell(i, j)` / `get_overlap(start, end)`: one scaled value.
- `get_values()`: the whole scaled array.
- `view(grid)`: a numpy **view** (no copy) of a sub-grid, built with
  `OverlapGrid.slices`.

### `AbstractScorer` (frozen ABC)

The half of a condition that builds the map. It holds only map parameters and is
hashable, so equal scorers share one map.

- `build_map(seq, grid)`: pre-fills the values with 1, calls `_fill`, forces cells that
  run past the sequence end to 1, and returns an `OverlapMap`.
- `_fill(seq, grid, values)` *(abstract)*: writes the scores and returns the multiplier.

### `Condition` (dataclass, `eq=False`)

A scorer plus thresholds. Fields: `scorer`, `min_val`, `max_val`, `opt_val`, `name`,
`map` (set by `bind`), `version`.

- `__post_init__`: validates the thresholds, then sets `version` on the instance. That
  marks the end of construction.
- `_check_thresholds(...)` *(static)*: at least one threshold must be set.
- `__setattr__`: routes assignments to threshold fields through `update`, and blocks
  swapping `scorer` once the condition is bound.
- `update(min_val, max_val, opt_val)`: changes thresholds atomically and bumps `version`
  only if a value actually changed. `Overlaps` uses `version` to detect stale candidates.
- `grid`: the bound map's grid.
- `label`: `name`, or the scorer's class name if there's no name.
- `is_constraint`: `min_val` or `max_val` is set.
- `is_objective`: `opt_val` is set.
- `bind(overlap_map)`: attaches the map once. Called by `Overlaps.add_condition`.
- `feasible()`: boolean mask of `min_val <= value <= max_val`. NaN counts as infeasible.
- `_plot_axes(ax)`: uses the given axes or creates new ones.
- `_plot_thresholds(ax)`: draws horizontal threshold lines.
- `plot_rows(a, b, ax)`: scatter of the minimum score per overlap start.

### `MapStore`

A cache of one map per scorer over one shared grid.

- `__init__(seq, grid)`: stores the sequence and grid, with
  `_maps: dict[AbstractScorer, OverlapMap]`.
- `get(scorer)`: builds the map on first use with `scorer.build_map` and makes it
  read-only.
- `view(scorer, grid)`: `get(scorer).view(grid)`. Called by `Overlaps.add_condition`.
- `__contains__`, `__len__`, `__repr__`: the usual container methods.

### `BackgroundKmers` (dataclass)

K-mer count table for a background. Fields: `kmer`, `keys` (sorted codes), `counts`.

- `__post_init__` / `_check_kmer(kmer)`: requires 2 ≤ kmer ≤ 16.
- `from_fasta(fasta, kmer)` / `from_seqs(seqs, kmer)` / `load(path)` *(classmethods)*:
  constructors.
- `save(path)`: writes an `.npz` file.
- `encode_seqs(seqs, kmer)` *(static)*: counts the k-mers on both strands, dropping and
  warning about invalid ones.
- `lookup(seq)`: the background count for each k-mer of `seq`, found with
  `searchsorted`.

### `KmerHomology(AbstractScorer)` (frozen, abstract)

Base for k-mer-based scorers. Field: `kmer=15`.

- `_kmer_scores(seq)` *(abstract)*: per-k-mer scores in 0..1, plus the multiplier.
- `_fill(...)`: worst k-mer score in each overlap, computed as a running max over
  sliding windows. Overlaps shorter than `kmer` score 0.

### `BGHomology(KmerHomology)`

Homology against a background genome. Fields: `cap`, `bg_fasta`, `bg_npz`, `bg`
(excluded from equality and the hash).

- `__post_init__`: checks `cap`, loads the background, checks that `kmer` matches, and
  lowers `cap` to the highest count or clips the counts to `cap`.
- `_load_bg()`: loads from FASTA if given, otherwise from the npz.
- `_kmer_scores(seq)`: `lookup(seq) / cap`, with multiplier `cap`.

### `SelfHomology(KmerHomology)`

Repeats within the sequence itself.

- `_kmer_scores(seq)`: builds `BackgroundKmers.from_seqs([seq])`, subtracts 1 from each
  count (the k-mer itself), and normalises by the maximum.

### `BoundaryMotif(AbstractScorer)`

Overlap opens and closes on a motif. Field: `motif="GC"`.

- `_fill(...)`: regex with a lookahead, so overlapping matches are all found. A cell
  scores 0 if a match starts at the overlap start and one ends at the overlap end,
  otherwise 1.

*(`PolyX` is commented out and marked TODO.)*

---

## `choppy/fragmentor.py`

### `NoAssemblyError(RuntimeError)`

Raised by `Fragmentor.chop` when no chain of overlaps spans the sequence.

### `GapError(NoAssemblyError)`

- `__init__(message, gaps)`: carries `gaps`, a list of `(start, end)` regions no
  segment can bridge. Raised by `Fragmentor.gap_check`.

### `Candidates`

Sorted parallel arrays of candidate overlaps: `starts`, `ends` and `extra_penalty`, all
read-only. They always include the sentinels `(0, 0)` and `(seq_len, seq_len)`.

- `__init__(seq_len, starts, ends, extra_penalty)`: sets up the sentinels, then calls
  `add`.
- `_freeze(arr)` *(static)*: makes an array read-only.
- `add(starts, ends, extra_penalty, skip_duplicates)`: validates, concatenates,
  lexsorts, and removes duplicate `(start, end)` pairs. It keeps the oldest one, or the
  newest when `skip_duplicates=False`.
- `__len__`: count including the sentinels.
- `n_overlaps`: count of real overlaps (without the sentinels).
- `__iter__`: `(start, end)` pairs as plain ints.
- `__repr__`.

### `Overlaps` (dataclass)

One level (segment or fragment) and its conditions and candidates. Fields: `seq`,
`min_len`, `max_len`, `min_step`, `conditions`, `candidates`, `maps` (a `MapStore`),
`grid` (an `OverlapGrid`), `manual` (a `Candidates` of hand-added overlaps), and
`_generated`, which holds the last generated candidates and the
`(condition, version)` state they were built from.

- `__post_init__`: builds its own `grid` and an empty `manual`. Uses the given `maps`
  (after checking that the sequence matches and the grid fits, via `grid.slices`) or
  creates its own `MapStore`. `min_step=1` becomes `None`.
- `shape`: the grid's shape.
- `constraints` / `objectives`: the conditions filtered by `is_constraint` /
  `is_objective`.
- `add_condition(condition)`: enforces a unique name (several conditions may share a scorer), then
  calls `condition.bind(maps.view(scorer, grid))`.
- `get_condition(name)`: looks up a condition by name.
- `feasible_mask()`: AND of `feasible()` over all constraints.
- `_thin(mask)`: drops the neighbours of each feasible cell within `min_step` (Manhattan
  distance in start/end).
- `generate_candidates(min_edge_len)`: feasible mask & `grid.edge_mask`, then `_thin`,
  then `grid.overlap`, then `Candidates`, then merges in `manual` (manual overlaps win).
  Records `_generated`.
- `_conditions_state()`: tuple of `(condition, version)` pairs, used for staleness.
- `candidates_stale`: true if there are no candidates, or if they were generated and
  the conditions have changed since. Candidates assigned directly are never stale.
- `add_candidates(starts, ends, extra_penalty)`: adds to `manual` and to the current
  `candidates`.
- `_require_candidates()`: returns `candidates` or raises.
- `_gaps()`: `(opens, closes)` per split. Each open is `start_i`; each close is a
  suffix-minimum of the ends.
- `gap_lengths()`: `closes - opens`.
- `plot_gap_lengths(min_gap_len, ax)`: histogram of the gap lengths.
- `plot_rows(indices, a, b)`: stacked subplots, one `Condition.plot_rows` per selected
  condition.
- `get_gapped_regions(min_gap_len)`: the `(open, close)` regions with length ≥ the
  threshold.

### `Assembly` (dataclass)

The result. Fields: `segments`, `overlaps`. `__len__` gives the number of segments.

### `Fragmentor` (dataclass)

The top-level API that owns both levels. Init params: `seq`, the segment overlap and
segment lengths, optional fragment overlap and fragment lengths, and `min_step`.
Fields: `grid` (the union of both levels' overlap ranges), `maps` (a shared
`MapStore`), `segment_overlaps` and `fragment_overlaps` (both `Overlaps`),
`use_fragments`, `seq_len`.

- `__post_init__(...)`: runs `check_seq`, validates both levels with `_check_lengths`,
  checks how fragments relate to segments, and builds the shared `grid`/`maps`. Then
  calls `_warn_unused_lengths` and builds the two `Overlaps` sharing `maps`.
- `_warn_unused_lengths(ranges)`: warns when most of the shared grid's lengths are used
  by neither level.
- `_check_lengths(name, ...)` *(static)*: checks types and positivity, that min ≤ max,
  and that max_len > min_len. Enforces `min_len > 2 * max_overlap_len`, lowering
  `max_overlap_len` with a warning or raising.
- `_level(fragment)`: returns the segment or fragment `Overlaps`, raising if there's no
  fragment level.
- `add_condition(condition, fragment)` / `get_condition(name, fragment)`: delegate to
  `_level`.
- `_lengths(fragment)`: `(min_len, max_len)` of the segments or fragments.
- `generate_candidates(fragment)`: `level.generate_candidates(min_len)`.
- `_current_candidates(fragment)`: regenerates first if `candidates_stale`.
- `plot_gap_lengths(...)` / `plot_rows(...)`: delegate to the level
  (`plot_gap_lengths` makes sure the candidates are current first).
- `gap_check(fragment, min_gap_len)`: uses `get_gapped_regions` with threshold
  `max_len + 1` (or `min_gap_len`) and raises `GapError` with the gaps.
- `add_candidates(..., fragment)`: delegates to the level.
- `chop()`: runs `gap_check` on both levels, then `_shortest_path` over the segment
  candidates. Raises `NoAssemblyError(_why_stuck(...))` if there's no path; otherwise
  `_walk_back` gives the `Assembly`.
- `_shortest_path(starts, ends, min_len, max_len)` *(static)*: one forward DAG sweep
  with unit edge weights. Each vertex relaxes a window found by `searchsorted`, keeping
  only lengths within the bounds. Returns `(weights, back)`.
- `_walk_back(back)` *(static)*: rebuilds the path from the back-pointers.
- `_why_stuck(starts, ends, weights)`: explains how far the chain got and how big the
  gap after that is.

---

## Major relations

```
Fragmentor
 ├─ grid: OverlapGrid (union of both levels' ranges)
 ├─ maps: MapStore(seq, grid) ──────────────── shared ─────────────┐
 ├─ segment_overlaps: Overlaps(maps=maps)                          │
 └─ fragment_overlaps: Overlaps(maps=maps) | None                  │
        ├─ grid: OverlapGrid (own range; must fit maps.grid.slices)│
        ├─ conditions: [Condition] ── bind(maps.view(scorer, grid))┘
        │        └─ scorer: AbstractScorer ── build_map → OverlapMap (cached in MapStore)
        │                 ├─ BoundaryMotif
        │                 └─ KmerHomology ─┬─ BGHomology ── BackgroundKmers (from fasta/npz)
        │                                  └─ SelfHomology ── BackgroundKmers.from_seqs
        ├─ manual: Candidates (hand-added, merged into every generation)
        └─ candidates: Candidates  ← generate_candidates (feasible_mask & edge_mask → _thin)
```

### Key flows

- **Map sharing:** each scorer's map is built once over the widest grid in `MapStore`.
  Each level reads it through a zero-copy `OverlapMap.view`, cut down with
  `OverlapGrid.slices`. Equal scorers get the same map because scorers are frozen and
  hashed by their parameters.
- **Staleness:** `Condition.update` and threshold assignments bump `version`.
  `Overlaps._conditions_state` compares the saved state with the current one, and
  `Fragmentor._current_candidates` regenerates when they differ.
- **Chop pipeline:** `chop` → `gap_check` (on both levels) → `_current_candidates` →
  `generate_candidates` → `_gaps` / `get_gapped_regions` → `_shortest_path` over the
  segment candidates → `_walk_back` → `Assembly`. Failures raise `GapError` or
  `NoAssemblyError`.
- **Fragments today:** the fragment level is generated and gap-checked, but `chop` only
  builds the path over the segment candidates. `objectives` and `extra_penalty` are
  stored but not yet used in the path (see the TODO in `_thin`).
