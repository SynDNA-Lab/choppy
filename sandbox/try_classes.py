# %%
from choppy.fragmentor import Fragmentor
from choppy.constraints import SelfHomology, BGHomology, BoundaryMotif
from Bio import SeqIO
import re
import numpy as np
# %%
all_sequences = list(SeqIO.parse("data/20240414-forSveta.fa", "fasta"))

seq_names = ["Maizel_COS-AT1G27430-GYF2", "Maizel_COS-SETH5", 
             "Pereira_COS-SynDNA-f1", "Pereira_COS-SynDNA-f2",
             "PV252688r_p6utr", "PQ537341r_p6utr,"
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
f.add_condition(SelfHomology(kmer=15, max_val=0.0))
f.add_condition(BGHomology(kmer=15, max_val=0.0, bg_fasta="data/S_cerevisiae-small.fa"))
f.add_condition(BoundaryMotif(motif="GC", min_val=0.5))
#asm = f.cut()   

# %%
f.segment_overlaps.constraints[2].plot_rows(0, 1000)
#%%
bm = BoundaryMotif(motif="GC", min_val=1.0)
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