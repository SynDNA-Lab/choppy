# %%
from choppy.fragmentor import Fragmentor
from choppy.constraints import SelfHomology, BGHomology, BoundaryMotif
from Bio import SeqIO

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
f.add_condition(BoundaryMotif(motif="GC", min_val=1.0))
asm = f.cut()   

#%%
bm = BoundaryMotif(motif="GC", min_val=1.0)