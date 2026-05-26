# %%
from Bio import SeqIO
from choppy.homology_finder import load_trie, create_kmer_trie, merge_tries, find_non_homologous_regions, find_local_non_homologous_regions, get_region_intersects
from choppy.primer_flanked_overlaps import get_primer_overlaps_from_seq, check_primer
import re
import primer3
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib.pyplot as plt

# %%
seqs = ["ACTGGGACACCAACAGTGTCCTCTAATTCTCACTCACATAGCGAGAAATCTAGTGTTTTCTACCAGCAA",
        "TTTTCTATCAACAGGAGTGGCCAGATAGTTATGCAACTGAAAAGGCTCTGAA",
        "AGTCTAGTTCTTACCCACAGAGGGAGAAGCCTAGTGTTTTGTACCCACAG",
        "ACCAGCTGACCAGATGACTGACACACCAGCAGTACCGTCTACTTTCTACTC",
        "AGCCCAGCATTTTCCACCAGCAGGCCTTGCCAGGTACTCATATACCTGAAGAGGC", 
        "TCCTGGACCAGTTGGCCAGACAACTGGCGCACCAACTATAACCTCTCCTTCCTACTCACAA"]

# %%

min_length = 17
max_length = 23
optimal_length = 20

# things that are easy to check: are presense of homopolymers and gc clamp
# so let's rely on that when scanning the regions

gc_clamp_pattern_forward = re.compile(r'[GC][AT][GC]|[AT][GC][GC]')
gc_clamp_pattern_reverse = re.compile(r'[GC][AT][GC]|[GC][GC][AT]')
poly_x_pattern = re.compile(r'(A{5,}|T{5,}|G{5,}|C{5,})')

# def check_primer(
#     primer, poly_x_pattern=None, poly_x_max_length=4,
#     min_gc=0.3, max_gc=0.7, 
#     min_tm=57, max_tm=62, 
#     max_hairpin_tm=24, max_homodimer_tm=45, 
#     max_3_self_tm=35.0
# ):
#     if poly_x_pattern is None:
#         poly_x_pattern = re.compile(r'(A{' + str(poly_x_max_length + 1) + r',}|T{' + str(poly_x_max_length + 1) + r',}|G{' + str(poly_x_max_length + 1) + r',}|C{' + str(poly_x_max_length + 1) + r',})')

#     if re.search(poly_x_pattern, primer):
#         return False, None, None
#     gc_content  = calc_gc_content(primer)
#     if gc_content < min_gc or gc_content > max_gc:
#         return False, gc_content, None
#     mt = primer3.calc_tm(primer)
#     if mt < min_tm or mt > max_tm:
#         return False, gc_content, mt
#     if primer3.calc_hairpin_tm(primer) > max_hairpin_tm:
#         return False, gc_content, mt
#     if primer3.calc_homodimer_tm(primer) > max_homodimer_tm:
#         return False, gc_content, mt
#     if primer3.calc_end_stability(primer, primer).tm > max_3_self_tm:
#         return False, gc_content, mt
#     return True, gc_content, mt

def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans('ATGC', 'TACG'))

calc_gc_content = lambda s: (s.count('G') + s.count('C')) / len(s)


# %%
for seq in seqs:
    for start in range(len(seq) - min_length + 1):
        for length in range(min_length, max_length + 1):
            if start + length > len(seq):
                continue
            primer = seq[start:start+length]
            is_valid, gc_content, tm = check_primer(primer, poly_x_max_length=4)
            if is_valid:
                print(f"Valid primer: {primer}, GC content: {gc_content:.2f}, Tm: {tm:.2f}")

# %%

filtered_primer_candidates = {}

for i in range(len(seqs)):
    seq_id = f"seq_{i}"
    print(seq_id)
    seq = str(seqs[i])
    cur_seq_forward_candidates = []
    cur_seq_reverse_candidates = []

    gc_matches = gc_clamp_pattern_forward.finditer(seq)

    for m in gc_matches:
        for pr_start in range(max(m.start() - (max_length - 3), 0), m.start() - min_length + 4 ):
            primer_cand = seq[pr_start:m.end()]
            valid, gc_content, mt = check_primer(primer_cand, poly_x_pattern)
            if not valid:
                continue
            cur_seq_forward_candidates.append({'seq': primer_cand, 'gc_content': gc_content, 'tm': mt, 'pos': (pr_start, m.end())})

    gc_matches = gc_clamp_pattern_reverse.finditer(seq)
    for m in gc_matches:
        for pr_end in range(m.end() + (min_length - 3), min(m.end() + (max_length - 3), len(seq)) + 1):
            primer_cand = reverse_complement(seq[m.start():pr_end])
            valid, gc_content, mt = check_primer(primer_cand, poly_x_pattern)
            if not valid:
                continue
            cur_seq_reverse_candidates.append({'seq': primer_cand, 'gc_content': gc_content, 'tm': mt, 'pos': (m.start(), pr_end)})

    filtered_primer_candidates[seq_id] = {"forward": cur_seq_forward_candidates, "reverse": cur_seq_reverse_candidates}

# %%
gc_clamp_pattern_relaxed = re.compile(r'[GC]')

filtered_primer_candidates_relaxed = {}

for i in range(len(seqs)):
    seq_id = f"seq_{i}"
    print(seq_id)
    seq = str(seqs[i])
    cur_seq_forward_candidates = []
    cur_seq_reverse_candidates = []

    gc_matches = gc_clamp_pattern_relaxed.finditer(seq)

    for m in gc_matches:
        for pr_start in range(max(m.start() - (max_length - 1), 0), m.start() - min_length + 2 ):
            primer_cand = seq[pr_start:m.end()]
            valid, gc_content, mt = check_primer(primer_cand, poly_x_pattern, verbose=True, max_tm=65.0)
            if not valid:
                continue
            cur_seq_forward_candidates.append({'seq': primer_cand, 'gc_content': gc_content, 'tm': mt, 'pos': (pr_start, m.end())})

    gc_matches = gc_clamp_pattern_relaxed.finditer(seq)
    for m in gc_matches:
        for pr_end in range(m.end() + (min_length - 1), min(m.end() + (max_length - 1), len(seq)) + 1):
            primer_cand = reverse_complement(seq[m.start():pr_end])
            valid, gc_content, mt = check_primer(primer_cand, poly_x_pattern, verbose=True, max_tm=65.0)
            if not valid:
                continue
            cur_seq_reverse_candidates.append({'seq': primer_cand, 'gc_content': gc_content, 'tm': mt, 'pos': (m.start(), pr_end)})

    filtered_primer_candidates_relaxed[seq_id] = {"forward": cur_seq_forward_candidates, "reverse": cur_seq_reverse_candidates}
