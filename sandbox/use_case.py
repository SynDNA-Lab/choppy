# %%
import warnings

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from choppy.constraints import (BackgroundKmers, BGHomology, BoundaryMotif,
                               Condition, SelfHomology)
from choppy.fragmentor import Fragmentor, GapError, NoAssemblyError

# %%
all_sequences = {rec.id: rec for rec in SeqIO.parse("data/20240414-forSveta.fa", "fasta")}

record = all_sequences["Pereira_COS-SynDNA-f2"]
seq = str(record.seq).upper()
print(record.id, len(seq), "bp")
# %%
f = Fragmentor(
    seq,
    min_segment_overlap_len=60, max_segment_overlap_len=120,
    min_segment_len=1500, max_segment_len=3000,
    min_fragment_overlap_len=30, max_fragment_overlap_len=60,
    min_fragment_len=500, max_fragment_len=1000,
    min_step=10,
)
# %% [markdown]
bg_kmers = BackgroundKmers.load("data/S_cerevisiae-small-k15.npz")

print(f"{len(bg_kmers.keys):,} distinct 15-mers in the background")

# %%
f.add_condition(Condition(SelfHomology(kmer=15), max_val=0.0, name="self"))
f.add_condition(Condition(BGHomology(kmer=15, cap=1024, bg_npz="data/S_cerevisiae-small-k15.npz"),
                          max_val=0.0, name="bg"))
f.add_condition(Condition(BoundaryMotif(motif="[GC]{3}"), max_val=0.0, name="gc3"))

f.segment_overlaps.conditions

print("constraints:", [c.label for c in f.segment_overlaps.constraints])
print("objectives: ", [c.label for c in f.segment_overlaps.objectives])
print("maps built: ", len(f.maps))

# %%
f.get_condition("bg").plot_rows()

# %%
f.plot_rows(a=0, b=4000)   
f.plot_rows(indices=["self", "gc3"], a=14000, b=19000) 

# %%
candidates = f.generate_candidates()
print(candidates, "->", candidates.n_overlaps, "real overlaps (plus the two sentinels)")
print("first few:", [ov for ov in list(candidates)[:5]])

# %%
f.plot_gap_lengths()
gaps = f.segment_overlaps.gap_lengths()
print("gaps: median", int(np.median(gaps)), "max", int(gaps.max()),
      "| max_segment_len", f.max_segment_len)

# %%
try:
    asm = f.chop()
    print(len(asm), "segments")
except GapError as e:
    print("GapError:", e)
    stuck = e.gaps
    (a, b), = stuck
except NoAssemblyError as e:
    print("NoAssemblyError:", e)

f.plot_rows(indices=["self", "gc3"], a=a-1000, b=b+1000)
# %%
self_hom = f.get_condition("self")
rid = self_hom.add_relaxed(max_val=4.0, extra_penalty=1.0, start=a, end=b)
print("relaxed version", rid, ":", self_hom.relaxed[-1])
print("candidates now stale:", f.segment_overlaps.candidates_stale)

# The plots read the condition live: the relaxed ceiling is drawn over the starts it
# applies to, labelled with its id, beside the main threshold.
f.plot_rows(indices=["self"], a=a-1000, b=b+1000)
# %%
asm = f.chop()
print(len(asm), "segments")
for i, (s, e) in enumerate(asm.segments, 1):
    print(f"  segment {i}: [{s:6d}, {e:6d})  {e - s:5d} bp")
print("junction overlaps:", asm.overlaps)

# %%
# Which junctions used the relaxed rule? Look up their penalty in the candidate list.
c = f.segment_overlaps.candidates
for s, e in asm.overlaps:
    i = int(np.flatnonzero((c.starts == s) & (c.ends == e))[0])
    flag = "  <- relaxed" if c.extra_penalty[i] else ""
    print(f"  [{s:6d}, {e:6d})  {e - s:3d} bp  penalty {c.extra_penalty[i]:.1f}{flag}")


# %%
# manually adding candidates
self_hom.remove_relaxed(rid)          # back to the strict rule, and a gap
f.generate_candidates()
print("gaps again:", f.segment_overlaps.get_gapped_regions(f.max_segment_len + 1))

f.add_candidates([16000], [16080], extra_penalty=0.5)
print("manual:", list(f.segment_overlaps.manual))
asm_manual = f.chop()
print(len(asm_manual), "segments with the hand-placed junction")

# %%
f.remove_candidates([16000], [16080])
print("manual now:", f.segment_overlaps.manual.n_overlaps,
      "| gaps:", f.segment_overlaps.get_gapped_regions(f.max_segment_len + 1))

# %% [markdown]
# edit the thresholds

print("before:", self_hom)
self_hom.update(max_val=2.0)          # tolerate 15-mers occurring up to twice more
f.get_condition("bg").max_val = 1.0   # and a little host-genome similarity
print("after: ", self_hom, "| stale:", f.segment_overlaps.candidates_stale)

asm_loose = f.chop()
print(len(asm_loose), "segments with the loosened rules:", asm_loose.segments)
# %%
# reusing the SelfHomology map that is already built (no new map, see `len(f.maps)`).
f.add_condition(Condition(SelfHomology(kmer=15), max_val=0.0, start=17000, end=len(seq),
                          name="self (strict)"))
print("maps built:", len(f.maps))
print("area rows:", f.get_condition("self (strict)").area_rows())

asm_area = f.chop()
print(len(asm_area), "segments:", asm_area.segments)

f.plot_rows(indices=["self", "self (strict)"]) 
# %% [markdown]
# min_step: how finely to sample the candidates

for step in (1, 10, 50):
    f.segment_overlaps.min_step = step
    print(f"min_step={step:3d}: {f.generate_candidates().n_overlaps:6d} candidates")

f.segment_overlaps.min_step = 10

# %%
f.explore(show=["self", "gc3", "bg"], a=13000, b=20000)