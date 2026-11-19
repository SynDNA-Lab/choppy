# %%
from choppy.fragmentor import Fragmentor, GapError, NoAssemblyError
from choppy.constraints import Condition, SelfHomology, BGHomology, BoundaryMotif
from Bio import SeqIO
import re
import numpy as np
# %%
all_sequences = list(SeqIO.parse("data/20240414-forSveta.fa", "fasta"))

seq_names = ["Maizel_COS-AT1G27430-GYF2", "Maizel_COS-SETH5",
             "Pereira_COS-SynDNA-f1", "Pereira_COS-SynDNA-f2",
             "PV252688r_p6utr", "PQ537341r_p6utr",
             "PQ488556r_p6utr", "PX021458r_p6utr",
             "MZ289137_rep", "OR500095r_naive"]

sequences = [seq for seq in all_sequences if seq.id in seq_names]

for seq in sequences:
    seq.seq = seq.seq.upper()

sequences_by_id = {seq.id: seq for seq in sequences}

seq = str(sequences[0].seq)
# %%
f = Fragmentor(seq,
               min_segment_overlap_len=50, max_segment_overlap_len=100,
               min_segment_len=200, max_segment_len=1000,
               min_step=10)

# A Condition = a scorer (builds the map, shared between levels) + thresholds.
# Thresholds are in the scorer's own units (repeat counts, background counts, 0/1 for motifs).
f.add_condition(Condition(SelfHomology(kmer=15), max_val=0.0, name="self"))
f.add_condition(Condition(BGHomology(kmer=15, cap=1024, bg_fasta="data/S_cerevisiae-small.fa"),
                          max_val=1.5, name="bg"))
f.add_condition(Condition(BoundaryMotif(motif="GC"), max_val=0.0, name="gc"))

f.segment_overlaps.constraints

# %%
# A single condition: by name, or by position in f.segment_overlaps.constraints
f.get_condition("bg").plot_rows()

# %%
# Several constraints stacked, sharing the x axis; names and indices can be mixed
f.plot_rows(a=0, b=1000)
f.plot_rows(indices=["self", 2], a=0, b=1000)

# %%
f.generate_candidates()
print(f.segment_overlaps.candidates)
f.plot_gap_lengths()

# %%
try:
    asm = f.chop()
    print(len(asm), "segments:", asm.segments)
except GapError as e:
    print(e)
except NoAssemblyError as e:
    print(e)

# %%
# Update thresholds: either update(...) (None removes one, omitted ones are kept)
# or plain assignment. Both bump the condition's version, which makes the candidates stale.
self_hom = f.get_condition("self")
print("before:", self_hom, "| stale:", f.segment_overlaps.candidates_stale)

self_hom.update(max_val=1.0)
f.get_condition("bg").max_val = 3.0
print("after: ", self_hom, "| stale:", f.segment_overlaps.candidates_stale)

# %%
# The plots read the thresholds live, so they show the new lines
f.plot_rows(indices=["self", "bg"], a=0, b=1000)

# %%
# gap_check / chop / plot_gap_lengths regenerate stale candidates on their own
f.plot_gap_lengths()
print(f.segment_overlaps.candidates, "| stale:", f.segment_overlaps.candidates_stale)

# %%
try:
    asm = f.chop()
    print(len(asm), "segments:", asm.segments)
except GapError as e:
    print(e)
except NoAssemblyError as e:
    print(e)

#%%
bm = BoundaryMotif(motif="GC", max_val=0.0)
bm.build_map(seq, min_len=50, max_len=100)

matches = list(re.finditer(bm.motif, seq))
starts = np.array([m.start() for m in matches], dtype=np.int64)
ends = np.array([m.end() for m in matches], dtype=np.int64)

opens = np.zeros(bm.map.shape[0], dtype=bool)
opens[starts[starts < bm.map.shape[0]]] = True

closes = np.zeros(len(seq) + 1, dtype=bool)
closes[ends] = True

i = np.arange(bm.map.shape[0])[:, None]
j = np.arange(bm.map.shape[1])[None, :]
min_len = 50
stop = i + min_len + j
shuts = closes[np.minimum(stop, len(seq))] & (stop <= len(seq))

# %%
opens[starts[starts < bm.map.shape[0]]]

# %%
min_step = 10

offsets = []
for i in range(min_step // 2):
    for j in range(min_step - i * 2):
        offsets.append((i, j))
offsets.pop(0)  # don't include (0, 0) because it would always be kept
offsets = np.array(offsets) 