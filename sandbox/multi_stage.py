# %%
from Bio import SeqIO
from choppy.homology_finder import load_trie, create_kmer_trie, merge_tries, find_non_homologous_regions, find_local_non_homologous_regions, get_region_intersects
from choppy.primer_flanked_overlaps import get_primer_overlaps_from_seq
import re
import primer3
import matplotlib.pyplot as plt
import networkx as nx
import heapq
from collections import defaultdict
import matplotlib.pyplot as plt
# %%
sequences = list(SeqIO.parse("data/random_seqs_50k.fa", "fasta"))

test_seq = str(sequences[0].seq)
kmer_size = 15
primer_kmer_size = 10
bg_trie = load_trie("../data/S_cerevisiae-R64-GCA_000146045_cat_15.marisa")
seq_tries = defaultdict(dict)

for seq in sequences:
    print(f"Processing {seq.id}...")
    trie = create_kmer_trie(seq, primer_kmer_size, bg=False)
    seq_tries[seq.id] = trie

full_seq_trie = merge_tries(list(seq_tries.values()))

if kmer_size != primer_kmer_size:
    for seq in sequences:
        print(f"Processing {seq.id}...")
        trie = create_kmer_trie(seq, kmer_size, bg=False)
        seq_tries[seq.id] = trie
    

# %%
threshold = 50
min_primer_length = 17
max_primer_length = 23
pr_threshold = min_primer_length
neighbourhood_size = 5500

bg_regions = {}
full_seq_regions = {}
cur_seq_regions = {}
local_neigh_regions = {}
for seq in sequences:
    bg_regions[seq.id] = find_non_homologous_regions(seq, bg_trie, [], kmer_size, threshold=threshold)
    full_seq_regions[seq.id] = find_non_homologous_regions(seq, full_seq_trie, [], kmer_size, threshold=pr_threshold)
    cur_seq_regions[seq.id] = find_non_homologous_regions(seq, seq_tries[seq.id], bg_trie, kmer_size, threshold=threshold)
    local_neigh_regions[seq.id] = find_local_non_homologous_regions(seq, kmer_size, threshold, neighbourhood_size)

# %%
primer_candidate_regions = {}
for seq in sequences:
    local_neigh_regions[seq.id] = get_region_intersects(local_neigh_regions[seq.id], bg_regions[seq.id], threshold)
    primer_candidate_regions[seq.id] = get_region_intersects(full_seq_regions[seq.id], local_neigh_regions[seq.id], pr_threshold)

# %%
primer_flanked_overlaps = {}

gc_clamp_pattern_forward = re.compile(r'[GC][AT][GC]|[AT][GC][GC]')
gc_clamp_pattern_reverse = re.compile(r'[GC][AT][GC]|[GC][GC][AT]')

for seq in sequences:
    print(seq.id)
    primer_flanked_overlaps[seq.id] = get_primer_overlaps_from_seq(str(seq.seq), primer_candidate_regions[seq.id], 
                                                                   gc_clamp_pattern_forward, gc_clamp_pattern_reverse, 3,
                                                                   five_prime_clamp_pattern=re.compile(r'[GC]'), min_overlap=50, 
                                                                   min_length=min_primer_length, max_length=max_primer_length)

# %%
min_step = 50
min_overlap = 50
max_overlap = 100

segment_overlaps = {}
for seq in sequences:
    print(seq.id)
    cur_overlaps = []
    for region in cur_seq_regions[seq.id]:
        if region[1] - region[0] < max_overlap and region[1] - region[0] >= min_overlap:
            cur_overlaps.append(region)
            continue
        for start in range(region[0], region[1] - max_overlap + 1, min_step):
            cur_overlaps.append((start, start + max_overlap))
        for end in range(region[1], region[0] + max_overlap - 1, -min_step):
            cur_overlaps.append((end - max_overlap, end))
    segment_overlaps[seq.id] = sorted(cur_overlaps, key=lambda x: x[0])

# %%
def get_edge_penalty(left_ov, right_ov, 
                     opt_primer_len, max_primer_len_diff, 
                     opt_overlap_len, max_overlap_len_diff, 
                     max_tm_diff=5.0):
    penalty = 0.0
    if left_ov['type'] == 'fr' and right_ov['type'] == 'fr':
        tm_diff = abs(left_ov['tm_left'] - right_ov['tm_right'])
        penalty += tm_diff / max_tm_diff * 0.1
    
    if left_ov['type'] == 'fr':
        primer_len_diff = abs(left_ov['pr_left_len'] - opt_primer_len)
        penalty += primer_len_diff / max_primer_len_diff * 0.1
    if right_ov['type'] == 'fr':
        primer_len_diff = abs(right_ov['pr_right_len'] - opt_primer_len)
        penalty += primer_len_diff / max_primer_len_diff * 0.1
    
    left_overlap_len = left_ov['pos'][1] - left_ov['pos'][0]
    right_overlap_len = right_ov['pos'][1] - right_ov['pos'][0]
    penalty += abs(left_overlap_len - opt_overlap_len) / max_overlap_len_diff * 0.1
    penalty += abs(right_overlap_len - opt_overlap_len) / max_overlap_len_diff * 0.1

    return round(penalty * 100)

def get_edge_clamp(left_ov, right_ov):
    return left_ov['left_clamp'] + right_ov['right_clamp']
    
def add_vertex(G, ov, overlap_list=None):
    if not G.has_node(ov['pos']):
        if ov['type'] == 'seg':
            G.add_node(ov['pos'], type=ov['type'], pos=ov['pos'],
                       left_clamp=0, right_clamp=0)
        if ov['type'] == 'fr':
            left_clamp = 3 if overlap_list[ov['pr_ind']]['forward']['seq'][0] == 'G' else 2
            right_clamp = 3 if overlap_list[ov['pr_ind']]['reverse']['seq'][0] == 'G' else 2
            G.add_node(ov['pos'], type=ov['type'], pos=ov['pos'],
                       pr_left_len=len(overlap_list[ov['pr_ind']]['forward']['seq']), 
                       pr_right_len=len(overlap_list[ov['pr_ind']]['reverse']['seq']), 
                       tm_left=overlap_list[ov['pr_ind']]['forward']['tm'], 
                       tm_right=overlap_list[ov['pr_ind']]['reverse']['tm'],
                       left_clamp=left_clamp, right_clamp=right_clamp)

# %%
seq_id = sequences[0].id

G = nx.DiGraph()

max_heterodimer_tm = 45.0
max_heterodimer_end_stability = 35.0

min_length = 300
max_length = 1800

opt_primer_len = (max_primer_length + min_primer_length) / 2
max_primer_len_diff = max_primer_length - min_primer_length

start_vertex = (-100, 0)
end_vertex = (len(str(sequences[0].seq)), len(str(sequences[0].seq)) + 100)

all_overlaps = (
    [{'pos': ov, 'type': 'seg'} for ov in segment_overlaps[seq_id]] +
    [{'pos': (ov['overlap_start'], ov['overlap_end']), 'type': 'fr', 'pr_ind': i} for i, ov in enumerate(primer_flanked_overlaps[seq_id])] +
    [{'pos': start_vertex, 'type': 'seg'}, 
     {'pos': end_vertex, 'type': 'seg'}]
)
all_overlaps = sorted(all_overlaps, key=lambda x: x['pos'][0])

#not_passed_count = 0
for i, ov in enumerate(all_overlaps):

    if i % 1000 == 0:
#        print(f"Processing overlap {i}/{len(all_overlaps)}, not_passed_count: {not_passed_count}")
        print(f"Processing overlap {i}/{len(all_overlaps)}")

    add_vertex(G, ov, primer_flanked_overlaps[seq_id])
    j = i + 1
    too_far = False
    allow_seg = ov['type'] == 'fr'
    while j < len(all_overlaps) and not too_far:
        next_ov = all_overlaps[j]
        if next_ov['pos'][0] - ov['pos'][0] > max_length:
            too_far = True
            continue
        length = next_ov['pos'][1] - ov['pos'][0]
        if min_length <= length <= max_length and (allow_seg or next_ov['type'] == 'fr'):
            add_vertex(G, next_ov, primer_flanked_overlaps[seq_id])
            #hetero_tm = 0.0
            #hetero_end_stability = 0.0
            # if ov['type'] == 'fr' and next_ov['type'] == 'fr':
                # they never got triggered anyway, but take too much time to compute for all pairs
                # hetero_tm = primer3.calc_heterodimer_tm(
                #     primer_flanked_overlaps[seq_id][ov['pr_ind']]['forward']['seq'], 
                #     primer_flanked_overlaps[seq_id][next_ov['pr_ind']]['reverse']['seq']
                # )
                # hetero_end_stability = primer3.calc_end_stability(
                #     primer_flanked_overlaps[seq_id][ov['pr_ind']]['forward']['seq'], 
                #     primer_flanked_overlaps[seq_id][next_ov['pr_ind']]['reverse']['seq']
                # ).tm
            G.add_edge(ov['pos'], next_ov['pos'], 
                       weight=100 + get_edge_penalty(G.nodes[ov['pos']], G.nodes[next_ov['pos']], 
                                                   opt_primer_len, max_primer_len_diff, 
                                                   max_overlap, max_overlap - min_overlap))
            #not_passed_count += 0 if G.edges[ov['pos'], next_ov['pos']]['hetero_pass'] else 1           
        j += 1

# %%
num_nodes = G.number_of_nodes()
print(f"Number of vertices: {num_nodes}")

num_edges = G.number_of_edges()
print(f"Number of edges: {num_edges}")

avg_in_degree = G.number_of_edges() / G.number_of_nodes()
avg_total_degree = (2 * G.number_of_edges()) / G.number_of_nodes()

print(f"Average In/Out Degree: {avg_in_degree}")
print(f"Average Total Degree: {avg_total_degree}")

# %%
edge_weights = [G.edges[e]['weight'] for e in G.edges()]
plt.figure(figsize=(10, 6))
plt.hist(edge_weights, bins=30, edgecolor='black')
plt.xlabel('Edge Weight')
plt.ylabel('Frequency')
plt.title('Histogram of Edge Weights')
plt.grid(True, alpha=0.3)
plt.show()

# %%
def reconstruct_path(predecessors, target_state):
    path = []
    step = target_state
    while step is not None:
        path.append(step[0]) # We only care about the vertex for the final path
        step = predecessors[step]
    path.reverse()
    return path


min_segment = 5000

in_state = (start_vertex, 0, 0) # (vertex, segment_length, clamp_code)
distances = defaultdict(lambda: float('inf'))
distances[in_state] = 0
predecessors = {in_state: None}

pq = []
heapq.heappush(pq, (0, in_state))
finished = False

rep_position = 0
rep_pos_step = 10000

while pq:
    cur_weight, cur_state = heapq.heappop(pq)
    if cur_weight > distances[cur_state]:
            continue
    
    if cur_state[0][0] > rep_position:
        print(f"At position {cur_state[0][0]}, weight: {cur_weight}")
        rep_position += rep_pos_step

    if cur_state[0] == end_vertex:
        print("Reached end vertex!")
        finished = True 
        break
    cur_vertex_data = G.nodes[cur_state[0]]

    for v in G.successors(cur_state[0]):
        edge_data = G.get_edge_data(cur_state[0], v)
        
        new_seg_len = round((cur_state[1] + v[1] - cur_state[0][1]) / 100) * 100
        new_clamp = cur_state[2]
        new_weight = cur_weight + edge_data['weight']

        v_data = G.nodes[v]

        # clamp-based restrictions
        if cur_state[2] == 0:
            # we have just started
            new_clamp = v_data['right_clamp']
        else:
            if new_seg_len > min_segment:
                # end of the segment, only left clamp of the current vertex matters (left of the new edge)
                if cur_vertex_data['left_clamp'] == 2 and new_clamp == 6: # C vs GG
                    continue
                elif cur_vertex_data['left_clamp'] == 3 and new_clamp == 4: # G vs CC
                    continue
                new_clamp = 0
            else:
                # in the middle of the segment, both clamps matter
                edge_clamp = cur_vertex_data['left_clamp'] + v_data['right_clamp']

                if cur_state[2] > 3:
                    if edge_clamp > 3 and cur_state[2] != edge_clamp:
                        continue
                    elif edge_clamp == 2 and cur_state[2] == 6: # C vs GG
                        continue
                    elif edge_clamp == 3 and cur_state[2] == 4: # G vs CC
                        continue
                else:
                    if cur_state[2] == 2 and edge_clamp == 6: # C vs GG
                        continue
                    elif cur_state[2] == 3 and edge_clamp == 4: # G vs CC
                        continue
        
                if new_clamp <= 3 and edge_clamp > 3:
                    new_clamp = edge_clamp
                #technically, this should never happen, since it requires a very small segment
                if new_clamp <= 3 and edge_clamp <= 3 and new_clamp != edge_clamp:
                    new_clamp += edge_clamp

        # distance-based restrictions
        if new_seg_len < min_segment:
            if v_data['type'] == 'seg':
                continue
        else:
            new_seg_len = v[1] - v[0]

        new_state = (v, new_seg_len, new_clamp)     

        if new_weight < distances[new_state]:
            distances[new_state] = new_weight
            predecessors[new_state] = cur_state
            heapq.heappush(pq, (new_weight, new_state))

# %%
min_length = 300 - 45
max_length = 1800 - 45
min_segment = 5500
segment_offset = 110
#offsets: for primers 45 in total, for segments 110 (from each side)

opt_primer_len = (max_primer_length + min_primer_length) / 2
max_primer_len_diff = max_primer_length - min_primer_length

start_vertex = (-100, 0)
end_vertex = (len(str(sequences[0].seq)), len(str(sequences[0].seq)) + 100)

all_overlaps = (
    [
        {'pos': ov, 'type': 'seg', 'left_clamp': 0, 'right_clamp': 0} 
        for ov in segment_overlaps[seq_id]
    ] +
    [
        {
            'pos': (ov['overlap_start'], ov['overlap_end']), 
            'type': 'fr', 
            'pr_ind': i,
            'left_clamp': 3 if ov['forward']['seq'][0] == 'G' else 2,
            'right_clamp': 3 if ov['reverse']['seq'][0] == 'G' else 2,
            'pr_left_len': len(ov['forward']['seq']),
            'pr_right_len': len(ov['reverse']['seq']),
            'tm_left': ov['forward']['tm'],
            'tm_right': ov['reverse']['tm']
        } 
        for i, ov in enumerate(primer_flanked_overlaps[seq_id])
    ] +
    [
        {'pos': start_vertex, 'type': 'seg', 'left_clamp': 0, 'right_clamp': 0}, 
        {'pos': end_vertex, 'type': 'seg', 'left_clamp': 0, 'right_clamp': 0}
    ]
)
all_overlaps = sorted(all_overlaps, key=lambda x: x['pos'][0])
overlap_states = [{} for _ in range(len(all_overlaps))]

overlap_states[0][(0, 0)] = (0, -1, None) # state: (segment_length, clamp_code): (weight, predecessor_ind, predecessor_state)

for i, ov in enumerate(all_overlaps):

    if i % 1000 == 0:
        print(f"Processing overlap {i}/{len(all_overlaps)}")

    j = i + 1
    too_far = False
    while j < len(all_overlaps) and not too_far:
        next_ov = all_overlaps[j]
        if next_ov['pos'][0] - ov['pos'][0] > max_length:
            too_far = True
            continue
        length = next_ov['pos'][1] - ov['pos'][0]
        if min_length <= length <= max_length:
            edge_weight = 100 + get_edge_penalty(ov, next_ov, 
                                                opt_primer_len, max_primer_len_diff, 
                                                max_overlap, max_overlap - min_overlap)
            for state, (cur_weight, _, _) in overlap_states[i].items():
                new_seg_len = round((state[0] + next_ov['pos'][1] - ov['pos'][1]) / 100) * 100
                new_clamp = state[1]
                new_weight = cur_weight + edge_weight

                if (new_seg_len >= min_segment or new_clamp == 0) and length + segment_offset > max_length:
                    continue
                if new_seg_len < min_segment and next_ov['type'] == 'seg':
                    continue

                # clamp-based restrictions
                if state[1] == 0:
                    # we have just started
                    new_clamp = next_ov['right_clamp']
                else:
                    if new_seg_len > min_segment:
                        # end of the segment, only left clamp of the current vertex matters (left of the new edge)
                        if ov['type'] == 'fr' and next_ov['type'] == 'fr':
                            if ov['right_clamp'] == 2 and next_ov['right_clamp'] == 6: # C vs GG
                                continue
                            elif ov['right_clamp'] == 3 and next_ov['right_clamp'] == 4: # G vs CC
                                continue
                        new_clamp = 0
                    else:
                        # in the middle of the segment, both clamps matter
                        edge_clamp = ov['right_clamp'] + next_ov['right_clamp']

                        if state[1] > 3:
                            if edge_clamp > 3 and state[1] != edge_clamp:
                                continue
                            elif edge_clamp == 2 and state[1] == 6: # C vs GG
                                continue
                            elif edge_clamp == 3 and state[1] == 4: # G vs CC
                                continue
                        else:
                            if state[1] == 2 and edge_clamp == 6: # C vs GG
                                continue
                            elif state[1] == 3 and edge_clamp == 4: # G vs CC
                                continue
            
                        if new_clamp <= 3 and edge_clamp > 3:
                            new_clamp = edge_clamp
                        #technically, this should never happen, since it requires a very small segment
                        if new_clamp <= 3 and edge_clamp <= 3 and new_clamp != edge_clamp:
                            new_clamp += edge_clamp

                if new_seg_len >= min_segment:
                    new_seg_len = 0
                
                new_state = (new_seg_len, new_clamp)
                if new_state not in overlap_states[j] or new_weight < overlap_states[j][new_state][0]:
                    overlap_states[j][new_state] = (new_weight, i, state)
        j += 1

            

        


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

# min_length = 18
# max_length = 23
# optimal_length = 20

# # things that are easy to check: are presense of homopolymers and gc clamp
# # so let's rely on that when scanning the regions

# gc_clamp_pattern_forward = re.compile(r'[GC][AT][GC]|[AT][GC][GC]')
# gc_clamp_pattern_reverse = re.compile(r'[GC][AT][GC]|[GC][GC][AT]')
# poly_x_pattern = re.compile(r'(A{5,}|T{5,}|G{5,}|C{5,})')

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

# def reverse_complement(seq):
#     return seq[::-1].translate(str.maketrans('ATGC', 'TACG'))

# calc_gc_content = lambda s: (s.count('G') + s.count('C')) / len(s)

# filtered_primer_candidates = {}

# for i in range(len(sequences)):
#     seq_id = sequences[i].id
#     print(seq_id)
#     seq = str(sequences[i].seq)
#     cur_seq_forward_candidates = []
#     cur_seq_reverse_candidates = []

#     for region in primer_candidate_regions[seq_id]:
#         gc_matches = gc_clamp_pattern_forward.finditer(seq, region[0], region[1])

#         for m in gc_matches:
#             for pr_start in range(max(m.start() - (max_length - 3), region[0]), m.start() - min_length + 4 ):
#                 primer_cand = seq[pr_start:m.end()]
#                 valid, gc_content, mt = check_primer(primer_cand, poly_x_pattern)
#                 if not valid:
#                     continue
#                 cur_seq_forward_candidates.append({'seq': primer_cand, 'gc_content': gc_content, 'tm': mt, 'pos': (pr_start, m.end())})

#         gc_matches = gc_clamp_pattern_reverse.finditer(seq, region[0], region[1])
#         for m in gc_matches:
#             for pr_end in range(m.end() + (min_length - 3), min(m.end() + (max_length - 3), region[1]) + 1):
#                 primer_cand = reverse_complement(seq[m.start():pr_end])
#                 valid, gc_content, mt = check_primer(primer_cand, poly_x_pattern)
#                 if not valid:
#                     continue
#                 cur_seq_reverse_candidates.append({'seq': primer_cand, 'gc_content': gc_content, 'tm': mt, 'pos': (m.start(), pr_end)})

#     filtered_primer_candidates[seq_id] = {"forward": cur_seq_forward_candidates, "reverse": cur_seq_reverse_candidates}

# # %%
# tms = [candidate['tm'] for candidate in filtered_primer_candidates['seq_0']['reverse']]
# plt.hist(tms, bins=20, edgecolor='black')
# plt.xlabel('Melting Temperature (°C)')
# plt.ylabel('Count')
# plt.title('Distribution of Filtered Primer Candidate Melting Temperatures')
# plt.show()

# # %%
# lengths = [len(candidate['seq']) for candidate in filtered_primer_candidates['seq_1']['forward']]
# plt.hist(lengths, bins=20, edgecolor='black')
# plt.xlabel('Primer Length (bp)')
# plt.ylabel('Count')
# plt.title('Distribution of Filtered Primer Candidate Primer Lengths')
# plt.show()

# # %%
# min_overlap = 60
# max_overlap = 100

# possible_overlaps = {}

# for seq_id in filtered_primer_candidates:
#     print(seq_id)
#     cur_possible_overlaps = []
#     for i, primer1 in enumerate(filtered_primer_candidates[seq_id]['forward']):
#         too_far = False
#         j = 0
#         while j < len(filtered_primer_candidates[seq_id]['reverse']) and not too_far:
#             primer2 = filtered_primer_candidates[seq_id]['reverse'][j]
#             if primer2['pos'][0] < primer1['pos'][1]:
#                 j += 1
#                 continue
#             if primer2['pos'][0] - primer1['pos'][1] > max_overlap:
#                 too_far = True
#                 continue
#             overlap_start = primer1['pos'][0]
#             overlap_end = primer2['pos'][1]
#             overlap_length = overlap_end - overlap_start
            
#             if min_overlap <= overlap_length <= max_overlap and primer3.calc_heterodimer_tm(primer1['seq'], primer2['seq']) <= 45:
#                 cur_possible_overlaps.append({
#                     'primer1': primer1,
#                     'primer2': primer2,
#                     'overlap_start': overlap_start,
#                     'overlap_end': overlap_end,
#                     'overlap_length': overlap_length
#                 })
#             j += 1
#     possible_overlaps[seq_id] = cur_possible_overlaps



# # %%
# primer_overlap_list = [(overlap['overlap_start'], overlap['overlap_end']) for overlap in possible_overlaps['seq_0']]

# lengths = [(len(prs["forward"]), len(prs["reverse"])) for prs in filtered_primer_candidates.values()]

# # %%
# list(filter(lambda x: x['pos'][0] == 414 or x['pos'][1] == 507, filtered_primer_candidates['seq_0']))

# # %%
# primer_cand = "AGGCTATGGGCCTCGCGA"

# rev_comp = primer_cand[::-1].translate(str.maketrans('ATGC', 'TACG'))

# primer3.calc_end_stability(primer_cand, primer_cand).tm

# # %%
# seq0_len = len(str(sequences[0].seq))
# covered = [False] * seq0_len

# for start, end in primer_overlap_list:
#     for i in range(start, end):
#         covered[i] = True

# gaps = []
# in_gap = False
# gap_start = 0

# for i in range(seq0_len):
#     if not covered[i] and not in_gap:
#         in_gap = True
#         gap_start = i
#     elif covered[i] and in_gap:
#         in_gap = False
#         gaps.append((gap_start, i, i - gap_start))

# if in_gap:
#     gaps.append((gap_start, seq0_len, seq0_len - gap_start))

# gap_lengths = [gap[2] for gap in gaps]
# plt.hist(gap_lengths, bins=20, edgecolor='black')
# plt.xlabel('Gap Length (bp)')
# plt.ylabel('Count')
# plt.title('Distribution of Uncovered Gap Lengths')
# plt.show()


# # %%
# primer_overlap_list = sorted([(overlap['overlap_start'], overlap['overlap_end']) for overlap in possible_overlaps['seq_0']], key=lambda x: x[0])

# # %%
# for seq in sequences:
#     seq_id = seq.id
#     primer_overlap_list = [(overlap['overlap_start'], overlap['overlap_end']) for overlap in possible_overlaps[seq_id]]

#     seq_len = len(str(seq.seq))
#     covered = [False] * seq_len

#     for start, end in primer_overlap_list:
#         for i in range(start, end):
#             covered[i] = True

#     gaps = []
#     in_gap = False
#     gap_start = 0

#     for i in range(seq_len):
#         if not covered[i] and not in_gap:
#             in_gap = True
#             gap_start = i
#         elif covered[i] and in_gap:
#             in_gap = False
#             gaps.append((gap_start, i, i - gap_start))

#     if in_gap:
#         gaps.append((gap_start, seq_len, seq_len - gap_start))

#     gap_lengths = [gap[2] for gap in gaps]
#     plt.hist(gap_lengths, bins=20, edgecolor='black')
#     plt.xlabel('Gap Length (bp)')
#     plt.ylabel('Count')
#     plt.title(f'Distribution of Uncovered Gap Lengths for {seq_id}')
#     plt.show()
