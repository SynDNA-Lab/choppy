# %%
from Bio import SeqIO
from choppy.homology_finder import load_trie, create_kmer_trie, merge_tries, find_non_homologous_regions, find_local_non_homologous_regions
import re
import primer3
import matplotlib.pyplot as plt

# %%
sequences = list(SeqIO.parse("data/random_seqs_50k.fa", "fasta"))

test_seq = str(sequences[0].seq)
kmer_size = 15
bg_trie = load_trie("../data/S_cerevisiae-R64-GCA_000146045_cat_15.marisa")
seq_tries = {}

for seq in sequences:
    print(f"Processing {seq.id}...")
    trie = create_kmer_trie(seq, kmer_size, bg=False)
    seq_tries[seq.id] = trie

full_seq_trie = merge_tries(list(seq_tries.values()))

# %%
threshold = 100
pr_threshold = 20
neighbourhood_size = 5500

bg_regions = {}
full_seq_regions = {}
cur_seq_regions = {}
local_neigh_regions = {}
for seq in sequences:
    bg_regions[seq.id] = find_non_homologous_regions(seq, bg_trie, [], kmer_size, threshold=threshold)
    full_seq_regions[seq.id] = find_non_homologous_regions(seq, full_seq_trie, bg_trie, kmer_size, threshold=pr_threshold)
    cur_seq_regions[seq.id] = find_non_homologous_regions(seq, seq_tries[seq.id], bg_trie, kmer_size, threshold=threshold)
    local_neigh_regions[seq.id] = find_local_non_homologous_regions(seq, kmer_size, threshold, neighbourhood_size)

# %%
def get_region_intersects(list1, list2, threshold):
    intersects = []
    
    j_start = 0
    for start1, end1 in list1:
        while j_start < len(list2) and list2[j_start][1] <= start1:
            j_start += 1
            
        for j in range(j_start, len(list2)):
            start2, end2 = list2[j]
            
            if start2 >= end1:
                break
            
            max_start = max(start1, start2)
            min_end = min(end1, end2)
            
            if max_start < min_end and (min_end - max_start) >= threshold:
                intersects.append((max_start, min_end))
                
    return intersects

local_bg_intersects = get_region_intersects(local_neigh_regions['seq_0'], bg_regions['seq_0'], threshold)
# primers should be inside valid overlap regions

primer_candidate_regions = {}
for seq in sequences:
    primer_candidate_regions[seq.id] = get_region_intersects(full_seq_regions[seq.id], local_neigh_regions[seq.id], pr_threshold)

# %%
# default primer 3 parameters
# length: 18-23, optimal: 20
# pimer tm: 57-62, optimal: 59 - max tm difference: 5
# gc content: 30-70%, optimal: 50%
# max 3' end stability: 9.0 (kcal/mol)
# max self complementarity: 45 (degrees?)
# max 3' end self complementarity: 35 (degrees?)
# max pair complementarity: 45 (degrees?)
# max 3' end pair complementarity: 35 (degrees?)
# max primer hairpin tm: 24 
# max poly-x: 4
# max gc in 3' end: 5

# monovalent cation concentration: 50 mM
# divalent cation concentration: 1.5 mM
# dNTP concentration: 0.6 mM
# anneling oligo concentration: 50.0 (something)

seq_id = 'seq_0'
min_length = 18
max_length = 23
optimal_length = 20

region = primer_candidate_regions[seq_id][0]

# things that are easy to check: are presense of homopolymers and gc clamp
# so let's rely on that when scanning the regions

gc_clamp_pattern = re.compile(r'[GC][AT][GC]|[AT][GC][GC]')
poly_x_pattern = re.compile(r'(A{5,}|T{5,}|G{5,}|C{5,})')


calc_gc_content = lambda s: (s.count('G') + s.count('C')) / len(s)

filtered_primer_candidates = []

for region in primer_candidate_regions[seq_id]:
    seq = str(sequences[0].seq)
    gc_matches = gc_clamp_pattern.finditer(seq, region[0] + min_length - 2, region[1])

    for m in gc_matches:
        for length in range(min_length, max_length + 1):
            pr_start = max(m.start() - (length - 3), region[0])
            primer_cand = seq[pr_start:m.end()]
            if re.match(poly_x_pattern, primer_cand):
                continue
            gc_content  = calc_gc_content(primer_cand)
            if gc_content < 0.3 or gc_content > 0.7:
                #print("GC content out of range for candidate:", primer_cand)
                continue
            mt = primer3.calc_tm(primer_cand)
            if mt < 57 or mt > 62:
                #print("Tm out of range for candidate:", primer_cand)
                continue
            if primer3.calc_hairpin_tm(primer_cand) > 24:
                #print("Hairpin tm too high for candidate:", primer_cand)
                continue
            if primer3.calc_homodimer_tm(primer_cand) > 45:
                #print("Homodimer tm too high for candidate:", primer_cand)
                continue
            filtered_primer_candidates.append({'seq': primer_cand, 'gc_content': gc_content, 'tm': mt, 'pos': (pr_start, m.end())})
 

# %%
tms = [candidate['tm'] for candidate in filtered_primer_candidates]
plt.hist(tms, bins=20, edgecolor='black')
plt.xlabel('Melting Temperature (°C)')
plt.ylabel('Count')
plt.title('Distribution of Filtered Primer Candidate Melting Temperatures')
plt.show()

# %%
min_overlap = 60
max_overlap = 100

possible_overlaps = []

for i, primer1 in enumerate(filtered_primer_candidates):
    too_far = False
    j = i + 1
    while j < len(filtered_primer_candidates) and not too_far:
        primer2 = filtered_primer_candidates[j]
        if primer2['pos'][0] - primer1['pos'][1] > max_overlap:
            too_far = True
            continue
        overlap_start = primer1['pos'][0]
        overlap_end = primer2['pos'][1]
        overlap_length = overlap_end - overlap_start
        
        if min_overlap <= overlap_length <= max_overlap and primer3.calc_heterodimer_tm(primer1['seq'], primer2['seq']) <= 45:
            possible_overlaps.append({
                'primer1': primer1,
                'primer2': primer2,
                'overlap_start': overlap_start,
                'overlap_end': overlap_end,
                'overlap_length': overlap_length
            })
        j += 1

# %%
overlap_list = [(overlap['overlap_start'], overlap['overlap_end']) for overlap in possible_overlaps]