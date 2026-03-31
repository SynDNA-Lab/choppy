# %%
from Bio import SeqIO
from choppy.homology_finder import load_trie, create_kmer_trie, merge_tries, find_non_homologous_regions, find_local_non_homologous_regions, get_region_intersects
from choppy.primer_flanked_overlaps import get_primer_overlaps_from_seq
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
local_bg_intersects = get_region_intersects(local_neigh_regions['seq_0'], bg_regions['seq_0'], threshold)
# primers should be inside valid overlap regions

primer_candidate_regions = {}
for seq in sequences:
    primer_candidate_regions[seq.id] = get_region_intersects(full_seq_regions[seq.id], local_neigh_regions[seq.id], pr_threshold)

# %%
primer_flanked_overlaps = {}

gc_clamp_pattern_forward = re.compile(r'[GC][AT][GC]|[AT][GC][GC]')
gc_clamp_pattern_reverse = re.compile(r'[GC][AT][GC]|[GC][GC][AT]')

for seq in sequences:
    print(seq.id)
    primer_flanked_overlaps[seq.id] = get_primer_overlaps_from_seq(str(seq.seq), primer_candidate_regions[seq.id], 
                                                                   gc_clamp_pattern_forward, gc_clamp_pattern_reverse, 3,
                                                                   five_prime_clamp_pattern=re.compile(r'[GC]'))


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

min_length = 18
max_length = 23
optimal_length = 20

# things that are easy to check: are presense of homopolymers and gc clamp
# so let's rely on that when scanning the regions

gc_clamp_pattern_forward = re.compile(r'[GC][AT][GC]|[AT][GC][GC]')
gc_clamp_pattern_reverse = re.compile(r'[GC][AT][GC]|[GC][GC][AT]')
poly_x_pattern = re.compile(r'(A{5,}|T{5,}|G{5,}|C{5,})')

def check_primer(
    primer, poly_x_pattern=None, poly_x_max_length=4,
    min_gc=0.3, max_gc=0.7, 
    min_tm=57, max_tm=62, 
    max_hairpin_tm=24, max_homodimer_tm=45, 
    max_3_self_tm=35.0
):
    if poly_x_pattern is None:
        poly_x_pattern = re.compile(r'(A{' + str(poly_x_max_length + 1) + r',}|T{' + str(poly_x_max_length + 1) + r',}|G{' + str(poly_x_max_length + 1) + r',}|C{' + str(poly_x_max_length + 1) + r',})')

    if re.search(poly_x_pattern, primer):
        return False, None, None
    gc_content  = calc_gc_content(primer)
    if gc_content < min_gc or gc_content > max_gc:
        return False, gc_content, None
    mt = primer3.calc_tm(primer)
    if mt < min_tm or mt > max_tm:
        return False, gc_content, mt
    if primer3.calc_hairpin_tm(primer) > max_hairpin_tm:
        return False, gc_content, mt
    if primer3.calc_homodimer_tm(primer) > max_homodimer_tm:
        return False, gc_content, mt
    if primer3.calc_end_stability(primer, primer).tm > max_3_self_tm:
        return False, gc_content, mt
    return True, gc_content, mt

def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans('ATGC', 'TACG'))

calc_gc_content = lambda s: (s.count('G') + s.count('C')) / len(s)

filtered_primer_candidates = {}

for i in range(len(sequences)):
    seq_id = sequences[i].id
    print(seq_id)
    seq = str(sequences[i].seq)
    cur_seq_forward_candidates = []
    cur_seq_reverse_candidates = []

    for region in primer_candidate_regions[seq_id]:
        gc_matches = gc_clamp_pattern_forward.finditer(seq, region[0], region[1])

        for m in gc_matches:
            for pr_start in range(max(m.start() - (max_length - 3), region[0]), m.start() - min_length + 4 ):
                primer_cand = seq[pr_start:m.end()]
                valid, gc_content, mt = check_primer(primer_cand, poly_x_pattern)
                if not valid:
                    continue
                cur_seq_forward_candidates.append({'seq': primer_cand, 'gc_content': gc_content, 'tm': mt, 'pos': (pr_start, m.end())})

        gc_matches = gc_clamp_pattern_reverse.finditer(seq, region[0], region[1])
        for m in gc_matches:
            for pr_end in range(m.end() + (min_length - 3), min(m.end() + (max_length - 3), region[1]) + 1):
                primer_cand = reverse_complement(seq[m.start():pr_end])
                valid, gc_content, mt = check_primer(primer_cand, poly_x_pattern)
                if not valid:
                    continue
                cur_seq_reverse_candidates.append({'seq': primer_cand, 'gc_content': gc_content, 'tm': mt, 'pos': (m.start(), pr_end)})

    filtered_primer_candidates[seq_id] = {"forward": cur_seq_forward_candidates, "reverse": cur_seq_reverse_candidates}

# %%
tms = [candidate['tm'] for candidate in filtered_primer_candidates['seq_0']['reverse']]
plt.hist(tms, bins=20, edgecolor='black')
plt.xlabel('Melting Temperature (°C)')
plt.ylabel('Count')
plt.title('Distribution of Filtered Primer Candidate Melting Temperatures')
plt.show()

# %%
lengths = [len(candidate['seq']) for candidate in filtered_primer_candidates['seq_1']['forward']]
plt.hist(lengths, bins=20, edgecolor='black')
plt.xlabel('Primer Length (bp)')
plt.ylabel('Count')
plt.title('Distribution of Filtered Primer Candidate Primer Lengths')
plt.show()

# %%
min_overlap = 60
max_overlap = 100

possible_overlaps = {}

for seq_id in filtered_primer_candidates:
    print(seq_id)
    cur_possible_overlaps = []
    for i, primer1 in enumerate(filtered_primer_candidates[seq_id]['forward']):
        too_far = False
        j = 0
        while j < len(filtered_primer_candidates[seq_id]['reverse']) and not too_far:
            primer2 = filtered_primer_candidates[seq_id]['reverse'][j]
            if primer2['pos'][0] < primer1['pos'][1]:
                j += 1
                continue
            if primer2['pos'][0] - primer1['pos'][1] > max_overlap:
                too_far = True
                continue
            overlap_start = primer1['pos'][0]
            overlap_end = primer2['pos'][1]
            overlap_length = overlap_end - overlap_start
            
            if min_overlap <= overlap_length <= max_overlap and primer3.calc_heterodimer_tm(primer1['seq'], primer2['seq']) <= 45:
                cur_possible_overlaps.append({
                    'primer1': primer1,
                    'primer2': primer2,
                    'overlap_start': overlap_start,
                    'overlap_end': overlap_end,
                    'overlap_length': overlap_length
                })
            j += 1
    possible_overlaps[seq_id] = cur_possible_overlaps



# %%
primer_overlap_list = [(overlap['overlap_start'], overlap['overlap_end']) for overlap in possible_overlaps['seq_0']]

lengths = [(len(prs["forward"]), len(prs["reverse"])) for prs in filtered_primer_candidates.values()]

# %%
list(filter(lambda x: x['pos'][0] == 414 or x['pos'][1] == 507, filtered_primer_candidates['seq_0']))

# %%
primer_cand = "AGGCTATGGGCCTCGCGA"

rev_comp = primer_cand[::-1].translate(str.maketrans('ATGC', 'TACG'))

primer3.calc_end_stability(primer_cand, primer_cand).tm

# %%
seq0_len = len(str(sequences[0].seq))
covered = [False] * seq0_len

for start, end in primer_overlap_list:
    for i in range(start, end):
        covered[i] = True

gaps = []
in_gap = False
gap_start = 0

for i in range(seq0_len):
    if not covered[i] and not in_gap:
        in_gap = True
        gap_start = i
    elif covered[i] and in_gap:
        in_gap = False
        gaps.append((gap_start, i, i - gap_start))

if in_gap:
    gaps.append((gap_start, seq0_len, seq0_len - gap_start))

gap_lengths = [gap[2] for gap in gaps]
plt.hist(gap_lengths, bins=20, edgecolor='black')
plt.xlabel('Gap Length (bp)')
plt.ylabel('Count')
plt.title('Distribution of Uncovered Gap Lengths')
plt.show()


# %%
primer_overlap_list = sorted([(overlap['overlap_start'], overlap['overlap_end']) for overlap in possible_overlaps['seq_0']], key=lambda x: x[0])

# %%
for seq in sequences:
    seq_id = seq.id
    primer_overlap_list = [(overlap['overlap_start'], overlap['overlap_end']) for overlap in possible_overlaps[seq_id]]

    seq_len = len(str(seq.seq))
    covered = [False] * seq_len

    for start, end in primer_overlap_list:
        for i in range(start, end):
            covered[i] = True

    gaps = []
    in_gap = False
    gap_start = 0

    for i in range(seq_len):
        if not covered[i] and not in_gap:
            in_gap = True
            gap_start = i
        elif covered[i] and in_gap:
            in_gap = False
            gaps.append((gap_start, i, i - gap_start))

    if in_gap:
        gaps.append((gap_start, seq_len, seq_len - gap_start))

    gap_lengths = [gap[2] for gap in gaps]
    plt.hist(gap_lengths, bins=20, edgecolor='black')
    plt.xlabel('Gap Length (bp)')
    plt.ylabel('Count')
    plt.title(f'Distribution of Uncovered Gap Lengths for {seq_id}')
    plt.show()
