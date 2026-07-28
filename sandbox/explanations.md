# Hierarchical Choppy

The algorithm is based on searching for a shortest path through a set of pre-generated overlaps, from a dummy start overlap `(0, 0)` to a dummy end overlap `(seq_len, seq_len)`. In a graph where overlaps are vertices, edges become fragments. Edges can be weighted to indicate "better" fragments. The shortest path will by design provide the smallest possible number of edges (with respect to weights).

*This logic is almost what choppy has now, but the current version weighs vertices (overlaps), not edges (fragments). Switching to edges gives easier control over fragment properties.*

## Overlap generation

For hierarchical mode, we construct two sets of overlaps: overlaps between segments, and overlaps between fragments within a segment. In non-hierarchical mode only segment overlaps should be used. The two may have independent sets of constraints, but the hard-coded difference is homology: fragment overlaps are checked for potential homology within their segment (+ background), segment overlaps are checked for homology against the entire sequence.

Some constraints act like hard thresholds, others add weight to an overlap or a fragment (any fragment using an overlap with `extra_weight > 0` will get its penalty added to its normal weight).

### Segment overlaps

The process is quite straightforward: we define suitable regions with choppy's `find_non_homologous_regions` and then get possible overlaps inside them.

### Fragment overlaps

Fragment overlaps are generated in several steps:

- For every position `i` of the sequence we define its range: how far we can go to the left or right without encountering the `seq[i:i+kmer_size]` kmer. (A fully unique kmer will then get the range `(0, seq_len)`.)

- Now, for any potential fragment overlap, we calculate its range by simply taking the intersection of the ranges of all its positions. At this stage only the unusably short overlap candidates are rejected:

    - `overlap[0] - range[0] < min_fragment_length` or `range[1] - overlap[1] < min_fragment_length` (there is no space to reach any valid segment overlap, or the overlap is outside of its own range).

    - `range[1] - range[0] < min_segment_length` (this overlap cannot be a part of any segment of valid length).

- In addition to that, fragment overlaps can be subject to any other constraints.

### TODO

- Define constraints as independent modules that can be applied to segment or fragment overlaps in any combination. *(A list of whatever we may need here will be helpful.)*
**Currently:** constraints are hardcoded.

- Design some filtering/grouping logic for almost identical overlaps. Generating too many almost identical candidates (like `(a, b)`, `(a + 1, b)`, `(a, b + 1)`, etc.) slows down the downstream algorithm significantly, increases memory usage and makes attempts to create a usable visualisation quite challenging.
**Currently:** there is some `min_step` logic (a new overlap must be at least `min_step` bp away unless something major changes, e.g. the range suddenly becomes much larger), but it is quite crude, and in some cases not very effective.
**Note:** we also need to take into account that sometimes a step back/forward can change things quite a bit — for instance, in the case of primer-flanked fragments, by resolving heterodimerisation or bringing the primer's melting temperature closer to (or into) a better range. There should be some trade-off between performance and the exhaustiveness of the search. Or some clever grouping logic that I haven't thought of yet.

- Decide on edge-case overlap constraint logic and other special cases. That's vague and is more related to the algorithm's to-do list. Stuff like restriction enzyme cut sites, or some cases where the normal constraints may be relaxed without extra costs... Anyway, let's say that this point is more about thinking through what can fall outside the normal scheme.

### Primer-flanked overlaps

That's quite a bit of code, so I put it here as a separate section, even though in the end it should be just a constraint module, applicable to both types of overlaps. There are lots of hardcoded things in here and some generalisation is clearly required.

- A full dictionary of potential 3' mispriming anchors is created (currently, all 6-mers are used).

- Primer search happens inside homology-free regions with respect to background (i.e. all potential fragment overlap locations).

- The search is based on going through all potential 3' ends of the primer (based on a GC-clamp regex). Every 3' end is checked against the mispriming anchors. If no matches are found, the primer is immediately considered misprime-free. If not, a melting temperature check is applied to all anchored positions. There are two misprime check modes:

    - **For pooled fragments:** mispriming is checked against any anchor in any of the sequences.

    - **For separated fragments:** mispriming is checked locally, `max_fragment_length` upstream or downstream of the primer. This ensures that the primer is unique within any fragment that its overlap can be a part of.

- The primer grows for the given 3' anchor, starting from `min_primer_length` and up to `max_primer_length`. Each primer is checked against multiple constraints. If it is clear that making the primer longer will make things worse (`max_tm` or `max_misprime_tm` exceeded), longer versions will not be checked.

Here is the block for primer constraints:
```python
# Check for single-nucleotide stretches (currently 5)
if cfg['poly_x_pattern'].search(primer_cand): break

# 5' clamp (was used as boundary motif but for fragments, currently any base)
if not cfg['five_prime_clamp_pattern'].search(primer_cand): continue

# GC-content (currently from 0.3 to 0.7)
gc_content = (primer_cand.count('G') + primer_cand.count('C')) / len(primer_cand)
if gc_content < cfg['min_gc'] or gc_content > cfg['max_gc']: continue

# Melting temperature (currently from 59 to 65)
tm = primer3.calc_tm(primer_cand, **cfg['tm_params'])
if tm > cfg['max_tm']: break
if tm < cfg['min_tm']: continue

# Melting temperature of a potential hairpin (currently 26 C)
if primer3.calc_hairpin_tm(primer_cand) > cfg['max_hairpin_tm']: continue
# Homodimer melting temperature (currently 47 C)
if primer3.calc_homodimer_tm(primer_cand) > cfg['max_homodimer_tm']: continue
# 3' homodimer melting temperature (currently 37 C)
if primer3.calc_end_stability(primer_cand, primer_cand).tm > cfg['max_3_self_tm']: continue

if side == "forward":
    native_3_prime_pos = m.end() - 1
    pos_tuple = (pr_start, m.end())
else:
    native_3_prime_pos = len(original_seq) - m.end()
    pos_tuple = (len(original_seq) - m.end(), len(original_seq) - pr_start)

# Mispriming melting temperature (currently 47 C)
if check_local_misprime(primer_cand, original_seq, native_3_prime_pos, side, seq_id, anchor_kmers, cfg):
    break
```

- After a set of forward and reverse primers is generated, overlaps are constructed by matching all allowed pairs, based on allowed overlap lengths and the range constraints for fragment overlaps described above.

#### TODO

The local mispriming constraint may be too strict. We do not allow mispriming within the longest possible fragment, but such an overlap could still be a part of a shorter fragment without any issues. This is similar to the logic of the homology-free ranges. However, it cannot be combined directly with ranges and, if implemented, will use its own check on the potential edges. It is easy to implement, but looks clumsy in cases when we don't use the primer-flanking constraint. So some generalised solution is required.

## Gap detection and constraint relaxation

After we have generated overlaps, we can already see some potentially problematic gaps (longer regions where there are no overlaps). If the gap is too long, it is clear that no path is possible — there is simply no fragment that can cover the entire gap. The potential solution is relaxed overlaps. We can add overlaps with some relaxed parameters and hard `extra_weight` penalties. Thus we ensure that they will be used either when absolutely necessary (`extra_weight` is prohibitively high — like 10–100 times the normal edge cost) or when their usage will help a lot (`extra_weight` is just high — like 1–4 times the normal edge cost).

### TODOs and challenges

- Currently, adding the relaxed parameters is done manually. There is a function to check for gaps, and a function to populate gaps with certain constraints relaxed, but the user must make this decision. I don't know how to automate it properly, since there is never a guarantee that a gap will be fixed and, therefore, there is a risk of an endless loop. Moreover, fully automated gap bridging will require a predefined logic of constraint importance. *We need to decide on the UI logic here.*

- Successfully bridging a gap does not guarantee successful path generation. Sometimes overlaps exist, but due to their short range some relaxed options are still required. This can be fixed with one of the following options (however, neither guarantees the result, and post-modifications may be required):

    - **Range-aware gap detection.** Currently, gaps are defined simply as regions without overlaps. This should be changed to more clever logic: first indicate overlaps that, due to their range, cannot reach the next one. That would be a recursive function, pruning the dead ends after the sets of overlaps have already been generated (this process must be aware of the current set of overlaps and cannot be based on the overlap range alone). The overlaps, however, should only be flagged, not removed, since they may become usable after the relaxed overlaps are added.

    **Pros:** useful in any case for downstream usability.
    **Cons:** can become complicated if we want to actually guarantee successful path generation. Also, it may still require multiple manually triggered iterations.

    - **Generous relaxed overlap seeding.** From the point of view of the shortest path algorithm, the relaxed overlaps will not be used unless needed, gap or no gap. Adding extra relaxed overlaps will not make the final result worse, only better.

    **Pros:** very straightforward logic from the user's point of view. They have to predefine a list of relaxed constraints and the corresponding `extra_weight`, and then everything else happens automatically.
    **Cons:** a computational nightmare. All relaxed overlaps will be considered by the algorithm and rejected only when a shorter, but otherwise equivalent, path is found. Potentially this will cause an explosion of states and will require some very clever logic to fix. Also, still no guaranteed solution.

## State generation and the shortest path

Once we are happy with the overlaps, we can move on to the shortest path search. We start from the leftmost overlap (`(0, 0)`) and move to the right. For each overlap we generate a dictionary of states, storing the lowest possible weight for the given state. A state is the information about the path taken to get to the given vertex that affects the choice of next steps. In non-hierarchical mode there is just one state: it doesn't matter how we got to this vertex, since the next steps are fully defined by the external parameters only. Thus a vertex stores only the lowest weight. In hierarchical mode, the state is a tuple of two values: `(segment_start, range_end)`.

- `segment_start` is the start position of the current segment. It defines whether we can (or have to) close the segment now, and how much extra penalty it will cost.

- `range_end` is the smallest among all `range[1]` values of all the fragment overlaps used so far in this segment. Since a segment must lie fully inside the intersection of the ranges of all its fragment overlaps, `range_end` defines the point before which we must close the segment.

The states are used as keys, and the corresponding value is `(weight, prev_vertex, prev_state)`:

- `weight` is the current total amount of accumulated penalties.

- `prev_vertex` is the index of the vertex from which we got here with this weight.

- `prev_state` is the state with which this step was made.

`weight` is being minimised; `prev_vertex` and `prev_state` don't affect the algorithm and are only stored to let us reconstruct the shortest path afterwards.

In the beginning, the first dummy vertex `(0, 0)` gets a single state `(0, seq_len)`.

When sitting at vertex `i`, the algorithm now checks all the available states one by one, and checks all possible next steps with the given state. For every vertex it can reach, it modifies the state dictionary, replacing the corresponding dictionary element if the current weight is smaller. Since we can only move from left to right, by the time the algorithm reaches vertex `i` it already stores a collection of all lowest weights for all the possible states. Therefore, after the last vertex `(seq_len, seq_len)` is reached, we pick the state with the lowest weight and go back, using `prev_vertex` and `prev_state`, thus obtaining the shortest path.

Basically, that's resource-constraint shortest path problem, but I stepped away from the implementations I could find.

### Notes, modifications and challenges

#### Weights

The weight that is going to be added at each step is calculated as a sum of the base edge weight (100), the calculated edge weight, the extra weight of the vertex (see the section on the relaxed overlaps) and the penalty for premature segment closing.

##### Edge penalty

The edge penalty is currently hardcoded and based on three parameters: the difference of primer melting temperatures (from 0 to 10), overlap length (from 0 to 10 for each overlap) and primer length (from 0 to 10 for each primer). Thus segment overlaps can only get an edge penalty from 0 to 20, based on the overlap lengths. Each of the penalties is based on optimal values, provided as external parameters, and changes linearly. `opt_primer_len` sets the optimal primer length (and with the current settings, extra long primers will actually get penalties higher than 10, since it is asymmetric); for the primer melting temperature difference the optimal value is 0, and for overlap length `max_overlap_length` is considered optimal.

##### Segment length

The segment length is regulated by `opt_segment_length` (currently 5000) and `min_segment_length` (currently 1000). `min_segment_length` is a hard threshold: the segment is not allowed to be closed if its length is lower. `opt_segment_length` acts as a non-strict threshold: when the segment length reaches it, the segment should be closed, but the resulting length at that point will be higher (for instance, if we are now at length 4900, the algorithm is allowed — and actually prompted — to add another full segment, thus potentially going up to ~5900: `opt_segment_length + max_fragment_length`). Premature segment closing is allowed as long as `min_segment_length` is reached, but adds an extra penalty `max(opt_seg_length - seg_len, 0) / opt_seg_length * 500`.

#### Offsets

In order to account for extra flanking sequences (like TT1 and TT2), offset parameters are introduced, with the current values:

```python
# Space for TT1
'segment_offset_left': 86,
'segment_offset_right': 116,
# Space for TT2
'seq_offset_left': 126,
'seq_offset_right': 105,
```

When the algorithm needs to start or close a segment, it subtracts the corresponding offset (or several offsets) from `max_fragment_length`. If we are starting or closing the sequence, the sequence offsets are also subtracted. In this way we leave enough space for the extra sequences, while at the same time allowing inner fragments to utilise the full available length range.

#### State dominance

**(!)** I've just realised that the idea is faulty and may lead to a missed solution, but I'll still describe it here.

In order to limit the amount of stored (and checked) states, and thus reduce runtime and memory usage, I've tried to use the idea that some states are better than others. If we reach a vertex with a lower weight and a larger `range_end` than an existing entry, we can safely replace the worse solution with the current one. A larger `range_end` includes all the possible steps of the lower one plus some more. I've applied the same logic to `seg_start`: as long as we're above `min_seg_length`, having a smaller segment length allows more options for the next steps (because every fragment overlap is checked against `seg_start`, since it has to lie within its range), and overshooting `opt_seg_length` doesn't cost extra penalty. But what I forgot is that premature closing will cost an extra penalty, which a segment with an earlier start may avoid. This approach may still be worth some thought and fixing, though.
However, using `range_end` alone reduces the number of stored states by only 1.2%, while the combination of `range_end` and `seg_start` reduced it by ~4 times.

It is probably fixable by taking into account the weight difference between the new state and the stored one. We know the size of the potential extra penalty due to premature closing. So if we accept only the weights that negate this extra penalty, it can work (unless I re-read it tomorrow and find a reason, why this is also problematic).

#### Performance

The main problem of the algorithm is the fast growing state space, since even a difference of 1 basepair generates an entirely independent path that has to be checked. On average each vertex gets ~80 states (early and bottleneck overlaps have fewer, overlaps that sit in the middle get considerably more). This requires some time and memory. The only way to fix it is to optimise the number of stored states (by careful rounding?) or the number of overlaps (by pre-choosing a better option among several similar ones).

I've already started that, but there may also still be some possibilities to optimise performance through careful caching.

### TODO

- Fix the weights. Currently some are hard coded, some are calculated based on the external parameters, some are defined by the external parameters. That's easier for me, since I have their scale of magnitude in my head, but is probably confusing for anyone else. And also just not great coding.

- Since I've added the `extra_weight` parameter for the overlap, it now makes sense to move the overlap length penalty from the edge to the vertex penalties.

- Refine the segment length penalty. The question here is how flexible we should be with overstepping `opt_segment_length`. Should it also be penalised, but allowed to grow much higher?

- Clearly separate fragment and segment overlaps, even if it leads to having duplicates. Currently the separation is somewhat smeared: a fragment overlap is allowed to act as a segment overlap if it satisfies all the required constraints. This is unnecessary and actually leads to some constraint hardcoding.

- Performance optimisation (?)

- To save memory (but not time), we may try to get rid of unused states after we've passed certain stages (like a segment end?) and can already say that some states will not be a part of the shortest path.

## Exploration and modifications

In the end, after the algorithm is done, we get a collection of states, among which there is a shortest path to each reached vertex (but not all possible valid paths). If the last vertex is reached, we can simply reconstruct the shortest path. But if it is not reached, things become more complicated. Of course it is easy to get the last reached vertex and the path to it, but to figure out why the path stopped I do manual checks, to see what the potential candidates in the region were, why they were rejected, and whether we can add something reasonable in the area. So far there is no system — and therefore no interface — for that, and this needs some work and consideration.

The things that can be done (or implemented) reasonably easily:

- get the shortest path to any reached vertex;

- recalculate the shortest path from any vertex (**(!)** it should either be forced to start a segment, or some way to provide — and explain — the start state is required). This will not fix the gap, but may provide two halves of a broken assembly;

- populate the hard region with some relaxed overlaps and rerun the algorithm.

## Interface

Overall, I imagine choppy (online or not) as a combination of two (optionally — one) sets of constraints that then define a hierarchical assembly. If it's online, it should be basically two sets of parameters instead of the current one. If it works as a package, then it is defined as a config or, probably, as a custom set of modules. This provides some (possibly incomplete) solution, after which manual changes may be required. So we also need functions for getting information on the current issues and fixing them. But here I still don't have a firm picture. This can be combined with some visual interface (even if it's not online), but I'm still unsure how it should look. An important thing is that the resulting states do not describe all the possibilities, only the optimal ones, so we can't just explore them. Some recalculations will be required.