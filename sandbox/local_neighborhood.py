# %%
from Bio import SeqIO
from collections import defaultdict 
import re

# %%
def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans('ATGC', 'TACG'))

# assumes that a is sorted
def next_val(a, val, default=None):
    return next((x for x in a if x > val), default)
def prev_val(a, val, default=None):
    return next((x for x in reversed(a) if x < val), default)

kmer_size = 15
optimal_overlap = 100
min_overlap = 50
min_fragment_size = 200
min_segment_size = 1000

min_step = 10

boundary_motif = re.compile(r'[GC]')
bm_len = 1
# %%
all_sequences = list(SeqIO.parse("data/20240414-forSveta.fa", "fasta"))

seq_id = "Maizel_COS-AT1G27430-GYF2"
seq = next(s for s in all_sequences if s.id == seq_id)

seq_str = str(seq.seq).upper()

# %%
kmers = defaultdict(list)
for i in range(len(seq_str) - kmer_size + 1):
    kmer = seq_str[i:i+kmer_size]
    kmers[kmer].append(i)
    kmers[reverse_complement(kmer)].append(i)

kmers = {k: v for k, v in kmers.items() if len(v) > 1}

# %%
homfree_ranges = [(0, len(seq_str))] * len(seq_str)

for i in range(len(seq_str)):
    start, end = homfree_ranges[i]
    if i < len(seq_str) - kmer_size + 1:
        kmer = seq_str[i:i+kmer_size]
        if kmer in kmers:
            n_val = next_val(kmers[kmer], i)
            if n_val is not None:
                end = min(end, n_val + kmer_size - 1)
    if i >= kmer_size - 1:
        kmer = seq_str[i-kmer_size+1:i+1]
        if kmer in kmers:
            p_val = prev_val(kmers[kmer], i - kmer_size + 1)
            if p_val is not None:
                start = max(start, p_val + 1)
    homfree_ranges[i] = (start, end)
    
# %%
def get_overlap_range(start, end, homfree_ranges, kmer_size):
    ov_end = min([el[1] for el in homfree_ranges[start:end - kmer_size + 1]])
    ov_start = max([el[0] for el in homfree_ranges[start + kmer_size - 1:end]])
    return ov_start, ov_end

def add_overlap(overlaps, cur_overlap):
    if cur_overlap["range"][0] != -1:
        overlaps.append(cur_overlap)

overlaps = []

prev_overlap = -min_step * 2
prev_range = (len(seq_str), 0)

for m in boundary_motif.finditer(seq_str):
    start_pos = m.span()[0]
    suboptimal_overlap = {"pos": (start_pos, start_pos + min_overlap), 
                          "range": (-1, len(seq_str) + 1) }
    for n in boundary_motif.finditer(seq_str, start_pos + min_overlap - bm_len):
        end_pos = n.span()[1]
        if end_pos - start_pos > optimal_overlap:
            add_overlap(overlaps, suboptimal_overlap)
            break
        cur_range = get_overlap_range(start_pos, end_pos, homfree_ranges, kmer_size)

        if prev_overlap != start_pos and start_pos - prev_overlap < min_step:
            if prev_range[0] <= cur_range[0] and prev_range[1] >= cur_range[1]:
                #print(f"Skipping candidate overlap: {start_pos}-{end_pos} because of previous one: {prev_overlap}")
                break
            
            prev_range = cur_range
            prev_overlap = start_pos
        if start_pos - prev_overlap >= min_step:
            prev_range = cur_range
            prev_overlap = start_pos

        if cur_range[1] - cur_range[0] < min_segment_size:
            add_overlap(overlaps, suboptimal_overlap)
            break
        
        if cur_range[0] > start_pos or cur_range[1] < end_pos:
            add_overlap(overlaps, suboptimal_overlap)
            break
        
        if cur_range[0] > suboptimal_overlap["range"][0] or cur_range[1] < suboptimal_overlap["range"][1]:
            add_overlap(overlaps, suboptimal_overlap)
        
        suboptimal_overlap = {"pos": (start_pos, end_pos), 
                              "range": cur_range}
    else:
        add_overlap(overlaps, suboptimal_overlap)

# %%
for ov in overlaps:
    reg = get_overlap_range(ov["pos"][0], ov["pos"][1], homfree_ranges, kmer_size)
    if reg != ov["range"]:
        print(f"Overlap {ov['pos']} has range {ov['range']} but should be {reg}")

# %%
list(filter(lambda reg: reg[0] > 0 or reg[1] < len(seq_str), homfree_ranges))