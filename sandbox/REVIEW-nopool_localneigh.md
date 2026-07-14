# Code review — `sandbox/nopool_localneigh.ipynb`

Scope: the primer-search cell, the homology-free-range helpers, and the shortest-path
cell. The size of the DP state space is a known issue and is deliberately left out.

Reviewed 2026-07-14, against the notebook as of commit `8aace73` (plus the config-dict
refactor of `get_overlap_states`). Findings are ordered by severity.

## Status

| # | Finding | Status |
|---|---|---|
| 1 | Reverse-primer coordinates in two frames | **fixed** (author), re-verified — one 1-base off-by-one left, see below |
| 2 | `check_misprime` dead, anchor map built over all sequences | **kept on purpose** — needed for future pooled designs |
| 3 | DP start sentinel can lose its sort position | **fixed** (author) — `sorted(key=lambda x: x['pos'])` |
| 4 | `get_edge_penalty` charges for a primer that will not exist | **partly fixed** (author) — residual, see below |
| 5 | Unasserted geometry invariants | **fixed** — `validate_config` + guard in `get_overlap_range` |
| 6 | `_has_boundary_motif` rescans the sequence per overlap | **fixed** — bounded match on both ends, using `seq_boundary_motif_len` |
| 7 | Magic numbers in the DP / penalty function | **deferred** |
| 8 | Minor / cleanup | **fixed** (except the `five_prime_clamp_pattern` no-op, kept as a config hook) |

Design intent confirmed with the author while reviewing:

- there is **no pooling** — each construct is amplified in its own reaction, so
  off-targets in other sequences do not matter;
- **segment ends carry no primers** — the boundary is handled by the TT1/TT2 vector
  arms (`segment_offset_left` / `segment_offset_right`);
- the optimal overlap is deliberately the longest allowed one;
- `opt_segment_length` is an "almost hard" maximum — overshooting by about one
  fragment is acceptable.

---

## 1. Reverse-primer coordinates live in two different frames (critical)

`build_anchor_kmers`, `find_primers_in_regions` and `check_local_misprime` disagree
about what `anchor_pos` / `native_3_prime_pos` mean on the reverse strand.

| function | value | frame |
|---|---|---|
| `find_primers_in_regions` | `native_3_prime_pos = len(original_seq) - m.end()` | 3' end, **original** coordinates |
| `build_anchor_kmers` | `len(seq_str) - pos - anchor_len` | anchor **start**, **revcomp** coordinates |
| `check_local_misprime` | compares the two, applies a distance window to them, and slices the revcomp string as if `anchor_pos` were a 3' end | mixed |

Forward primers are self-consistent; reverse primers are not. Measured on a 3 kb random
sequence: the on-target self-site is correctly excluded for **200/200 forward** primers
and **0/200 reverse** primers.

The damage is not limited to the missed self-exclusion. Because the mirrored coordinate
makes a reverse primer collide with *itself* whenever `L - native - k` happens to fall
inside `[native, native + max_frag_length]`, there is a systematic dead band of width
`max_frag_length / 2`, centred on the middle of every sequence, in which reverse primers
are discarded as misprimers. Measured on `Maizel_COS-AT1G27430-GYF2` (5925 nt, predicted
band 2459–2959):

| | reverse primers found | with 3' end inside the band |
|---|---|---|
| current code | 1618 | **2** |
| coordinates fixed | 1705 | **143** |

Every sequence in `data/20240414-forSveta.fa` has such a band (e.g. 12149–12649 for
`Pereira_COS-SynDNA-f2`). Outside the band the screen inspects off-target sites at
mirrored positions, so genuine near-downstream misprimers are never Tm-checked — the
fixed version rejects ~54 primers per sequence that the current code accepts.

**Fixed by the author**, by putting everything in the *original* frame instead: reverse
anchors are stored at `pos` (the 3' end of a reverse primer carrying that anchor),
`check_local_misprime` slices the original sequence, and the caller passes `original_seq`
rather than `search_seq`. Re-verified: all 11 850 entries of the anchor map now satisfy
the original-coordinate 3'-end convention, the self-site is excluded on both strands, and
the dead band is gone (143 reverse primers in the old dead zone of
`Maizel_COS-AT1G27430-GYF2`, versus 2 before).

The same edit also added the missing **reverse-complement misprime check** — off-target
sites are now screened on both strands, with the correct duplex partner in each case, and
the opposite-strand site is correctly required to lie *downstream* of the primer (where a
convergent spurious product could form).

**Still outstanding (1 base):** the forward fragment window is off by one relative to the
reverse one. `frag_start_pos = native_3_prime_pos - len(primer_cand)` puts the window
start one base 5' of the primer, whose 5' end is at `native - len + 1`; the reverse branch
computes its bound exactly. It should be `native_3_prime_pos - len(primer_cand) + 1`.

## 2. `check_misprime` is dead code — and it costs time

**Kept on purpose** — the cross-sequence screen and the all-sequence anchor map are wanted
for future pooled designs. The rest of this section is retained as documentation of what
it currently costs.

`check_misprime` is defined but never called; `find_primers_in_regions` only calls
`check_local_misprime`. Consequently `sequences_by_id` is threaded through
`find_primer_flanked_overlaps` → `find_primers_in_regions` → nowhere, and
`build_anchor_kmers(reverse_position=False)` is never used.

Given the no-pool design this is correct behaviour, so the function and its plumbing
should be deleted. There is a knock-on cost worth fixing at the same time:
`anchor_kmers` is built over **all 26 sequences**, and `check_local_misprime` then walks
each anchor's hit list discarding everything with `anchor_seq_id != seq_id`. About 25/26
of the entries it iterates are guaranteed misses. Building the anchor map per sequence
removes that work entirely.

## 3. The DP's start sentinel can lose its sort position (latent)

`get_all_overlaps` sorts by `x['pos'][0]` only. The start dummy is `(0, 0)`, and the sort
is stable with segment overlaps listed first — so any real overlap with `pos[0] == 0`
ties with the dummy and sorts *before* it. `get_overlap_states` hard-codes
`overlap_states[0][(0, seq_len)]`, so the DP would seed the wrong node and the true start
would be unreachable.

Not currently triggered, but only just: the earliest segment-overlap start is 1 in
`Maizel_COS-AT1G27430-GYF2` and 3 in `Pereira_COS-SynDNA-f2`. It fires as soon as a
sequence begins with `GC` inside a non-homologous region.

**Fixed by the author:** `sorted(all_overlaps, key=lambda x: x['pos'])`, which orders the
`(0, 0)` sentinel ahead of any real overlap starting at 0.

## 4. `get_edge_penalty` charges for a primer that will not exist

The edge weight is computed once per `(ov, next_ov)` pair, *before* the DP chooses
between the "start a new segment" and "continue the segment" branch. When `next_ov` is an
`fr` overlap that starts a new segment (the `boundary_motif` case), the fragment ends at
a segment boundary — where, per the design, the reverse primer is supplied by the TT arm,
not by the overlap. Yet the weight still includes `right_ov['pr_right_len']` and the
`tm_left`/`tm_right` difference against that non-existent primer.

**Partly fixed by the author:** `get_edge_penalty` gained a `primer_check` flag, passed
`False` on the seg-start edge and `True` on the continue edge, and the weight is now
computed inside the state loop. That correctly stops charging for the right-hand primer at
a segment boundary, and it makes the absent `is_heterodimer_risk` check in that branch
consistent rather than an oversight.

**Residual:** `primer_check=False` also drops the *left*-hand primer term, and the left
primer usually does exist. Which primers a fragment `(ov → next_ov)` really has:

- forward primer = `ov`'s — exists iff `ov['type'] == 'fr'` **and the fragment does not
  start the segment** (`state[0] != ov['pos'][0]`; if it does, the TT1 arm supplies it,
  which is exactly the case the DP already compensates for with `segment_offset_left`);
- reverse primer = `next_ov`'s — exists iff `next_ov['type'] == 'fr'` and the fragment does
  not end the segment (i.e. the continue branch).

So the flag should be two flags, and the left one depends on `state`, not just on the
branch. That is cheap now that the weight is computed inside the state loop. The
`tm_left`/`tm_right` difference term should only be charged when *both* primers exist.

## 5. Two unasserted invariants in the overlap/primer geometry

`get_overlap_range` slices `homfree_ranges[start:end - kmer_size + 1]` and calls `min()`
on it — an overlap shorter than `kmer_size` yields an empty slice and a `ValueError`.
Safe today (`min_overlap` 50, relaxed 40, vs `kmer_size` 15), but silently dependent on it.

`find_primer_flanked_overlaps` keeps a forward/reverse pair when
`rev.pos[1] - fwd.pos[0]` is in `[min_overlap, max_overlap]`. That both primers actually
lie *inside* the overlap they define holds only because
`min_overlap (50) > max_primer_length (30)`. Lower `min_overlap` below
`max_primer_length` and the code silently starts emitting overlaps with a primer hanging
outside the homology arm.

**Fixed:** added `validate_config(cfg)`, called at the top of
`find_primer_flanked_overlaps` so it covers the relaxed configs built by
`get_relaxed_overlaps` too. It rejects `min_overlap < max_primer_length`,
`min_overlap < kmer_size`, `min_overlap > max_overlap` and
`min_primer_length > max_primer_length`, each with a message saying what would break.
`get_overlap_range` additionally raises on an overlap shorter than `kmer_size` instead of
dying inside `min()`. Both the current config and the relaxed one (`min_overlap=40`) pass.

## 6. `_has_boundary_motif` rescans the sequence for every overlap

```python
ends_with_motif = any(m.end() == end for m in boundary_motif_pattern.finditer(seq, 0, end))
```

This scans from position 0 for every overlap — with ~10^4–10^5 overlaps on a 25 kb
sequence, ~10^9 scanned characters, when only the last few bases matter. `finditer` also
yields non-overlapping matches, so a motif ending exactly at `end` can be missed if it
overlaps an earlier match: harmless for the current 2 nt `GC` pattern, a real bug for any
self-overlapping motif (e.g. if the 3 nt clamp-style pattern were ever used here).

**Fixed:** `_has_boundary_motif` now takes `motif_len` (the new `seq_boundary_motif_len`
config key) and matches both ends in place — `pattern.match(seq, start, start + motif_len)`
and `pattern.match(seq, end - motif_len, end)`. `re.match` anchors at the given start and
the `endpos` keeps it inside the motif window, so a motif ending exactly at `end` can no
longer be shadowed by an earlier overlapping match. Verified identical to the old version
on 12 000 real overlaps across four sequences, and ~6x faster on that sample (the gap
widens with sequence length, since the old cost grew with `end`).

`get_all_overlaps` now takes `cfg` instead of the bare pattern, and reads both the pattern
and the length from it.

## 7. Magic numbers that belong in `CONFIG`

**Deferred by the author.**


Same family as the `opt_primer_len` / `max_primer_len_diff` TODO that was just closed:

- the base edge cost `100` and the short-segment penalty scale `500` in `get_overlap_states`;
- the four `0.1` term weights in `get_edge_penalty`;
- the call site passes `max_overlap` as `opt_overlap_len`, so "optimal overlap = longest
  allowed" is currently an accident of the argument list rather than a stated choice. The
  behaviour is intended — it just needs its own key.

Related, and **not** a bug: `opt_seg_length` acts as a hard cap in the continue branch
while the seg-start branch may overshoot it by one fragment. That matches the intended
"almost hard maximum", so it needs a comment, not a `max_segment_length`.

## 8. Minor / cleanup

**Fixed**, except where noted:

- `next_val` / `prev_val` were linear scans of the occurrence list, called twice per base in
  `compute_homfree_ranges`. The lists are already sorted, so they now use `bisect`
  (verified to agree with the old versions on 20 000 random cases).
- `reverse_complement` rebuilt its `maketrans` table on every call — including twice per
  base inside `build_anchor_kmers`. The table is now a module-level constant.
- `too_far = True; continue` in the DP is now a plain `break` (the author had already
  changed the `continue`); the redundant `too_far` flag and the `max(ov['pos'][0], 0)`
  guards are gone.
- `enumerate` in `get_all_overlaps` is dropped.
- `seg`-type overlaps and both sentinels now carry `boundary_motif: True` explicitly (true
  by construction for segment overlaps), so the DP condition no longer depends on
  short-circuit evaluation to avoid a `KeyError`.
- The "annotate with `range`, then filter by `min_frag_length`" block, copy-pasted three
  times, is now `annotate_overlap_ranges(overlaps, homfree_ranges, cfg)` (verified to
  reproduce the old block's output exactly).
- `_heterodimer_tm`'s `lru_cache` is bounded at 2^20 entries.
- **Not changed:** `five_prime_clamp_pattern` is `^[ATGC]`, so the check is currently a
  no-op, but it is a deliberate config hook — left in place.
- `_heterodimer_tm`'s `lru_cache` is unbounded.

## Checked and found sound

- The `break` on a positive local misprime (which discards all longer primers sharing that
  3' end): a longer primer contains every substring of the shorter one, so its best-alignment
  Tm against the same off-target window can only be higher. The pruning is valid.
- The `± (kmer_size - 1)` arithmetic in `compute_homfree_ranges` / `get_overlap_range`:
  the start and end bounds are consistent with "the duplicate k-mer must not be fully
  contained in the segment".
- Occurrence lists in `compute_homfree_ranges` are built in increasing order, so the
  `next_val` / `prev_val` sorted-input assumption holds. Palindromic k-mers self-collide
  into a `[i, i]` entry, but neither helper can return `i` itself, so it is harmless.

## Left to do

1. The 1-base off-by-one in the forward fragment window (finding 1).
2. The residual left-primer term in `get_edge_penalty` (finding 4).
3. Config keys for the DP's magic numbers (finding 7).

Of the changes made so far, only finding 1 moves the primer set — everything else is
behaviour-preserving on the current config (checked), or, in the case of finding 4, changes
edge weights. The notebook has not been re-run.
