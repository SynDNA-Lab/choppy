# %%
from Bio import SeqIO
from choppy.homology_finder import load_trie, create_kmer_trie, merge_tries, find_non_homologous_regions, find_local_non_homologous_regions


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

