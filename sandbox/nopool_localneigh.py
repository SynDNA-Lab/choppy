# %%
from Bio import SeqIO
from choppy.homology_finder import load_trie, create_kmer_trie, find_non_homologous_regions, find_local_non_homologous_regions, get_region_intersects
import primer3
from collections import defaultdict
import re
# %%
all_sequences = list(SeqIO.parse("data/20240414-forSveta.fa", "fasta"))

# seq_names = ["Maizel_COS-AT1G27430-GYF2", "Maizel_COS-SETH5", 
#              "Pereira_COS-SynDNA-f1", "Pereira_COS-SynDNA-f2",
#              "PV252688r_p6utr", "PQ537341r_p6utr,"
#              "PQ488556r_p6utr", "PX021458r_p6utr",
#              "MZ289137_rep", "OR500095r_naive"]

# sequences = [seq for seq in all_sequences if seq.id in seq_names]

sequences = all_sequences

for seq in sequences:
    seq.seq = seq.seq.upper()

sequences_by_id = {seq.id: seq for seq in sequences}

# %%
CONFIG = {
    'kmer_size': 15,
    'max_frag_length': 1000,
    'min_frag_length': 200,
    'opt_segment_length': 5000,
    'min_segment_length': 1000,
    'segment_offset_left': 78,
    'segment_offset_right': 108,
    'seq_offset_left': 106,
    'seq_offset_right': 105,
    'min_primer_length': 17,
    'max_primer_length': 30,
    'primer_3prime_anchor': 6,
    'min_step': 10,
    'min_gc': 0.3,
    'max_gc': 0.7,
    'min_tm': 57.0,
    'max_tm': 62.0,
    'tm_params': dict(
        mv_conc=222.0, dv_conc=0.0, dntp_conc=0.0, dna_conc=500.0,
        tm_method='santalucia', salt_corrections_method='schildkraut',
    ),
    'max_hairpin_tm': 24.0,
    'max_homodimer_tm': 45.0,
    'max_3_self_tm': 35.0,
    'max_misprime_tm': 47.0,
    'min_overlap': 50,
    'max_overlap': 100,
    'poly_x_pattern': re.compile(r'(A{5,}|T{5,}|G{5,}|C{5,})'),
    'clamp_length': 3,
    'clamp_pattern': re.compile(r'[GC][AT][GC]|[AT][GC][GC]'),
    'five_prime_clamp_pattern': re.compile(r'^[ATGC]')
}
# %%
bg_trie = load_trie("../data/S_cerevisiae-R64-GCA_000146045_cat_15.marisa")
seq_tries = {}
for seq in sequences:
    trie = create_kmer_trie(seq, CONFIG['kmer_size'], bg=False)
    seq_tries[seq.id] = trie 

bg_regions = {}
cur_seq_regions = {}
for seq in sequences:
    bg_regions[seq.id] = find_non_homologous_regions(seq, bg_trie, [], CONFIG['kmer_size'], threshold=CONFIG['min_overlap'])
    cur_seq_regions[seq.id] = find_non_homologous_regions(seq, seq_tries[seq.id], bg_trie, CONFIG['kmer_size'], threshold=CONFIG['min_overlap'])

# %%
def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans('ATGC', 'TACG'))

# helper search utilities (assume lists are sorted ascending)
def next_val(a, val, default=None):
    return next((x for x in a if x > val), default)

def prev_val(a, val, default=None):
    return next((x for x in reversed(a) if x < val), default)

def compute_homfree_ranges(seq_str, kmer_size):
    """Return a list of (start,end) homfree ranges for each base in seq_str.

    A position i is annotated with the nearest previous/next k-mer collision
    adjusted to k-mer coordinates, matching the original inline logic.
    """
    kmers = defaultdict(list)
    for i in range(len(seq_str) - kmer_size + 1):
        kmer = seq_str[i:i+kmer_size]
        kmers[kmer].append(i)
        kmers[reverse_complement(kmer)].append(i)

    kmers = {k: v for k, v in kmers.items() if len(v) > 1}

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

    return homfree_ranges

def build_anchor_kmers(sequences, anchor_len, reverse_position = True):
    """Builds the 3' k-mer hash map for fast off-target screening."""
    anchor_kmers = defaultdict(list)
    for seq in sequences:
        seq_str = str(seq.seq).upper()
        for pos in range(len(seq_str) - anchor_len + 1):
            kmer_fwd = seq_str[pos:pos + anchor_len]
            kmer_rev = reverse_complement(kmer_fwd)
            
            anchor_kmers[kmer_fwd].append((seq.id, pos + anchor_len - 1, "forward"))
            if reverse_position:
                anchor_kmers[kmer_rev].append((seq.id, len(seq_str) - pos - anchor_len, "reverse"))
            else:
                anchor_kmers[kmer_rev].append((seq.id, pos + anchor_len - 1, "reverse"))
            
    return anchor_kmers

def check_misprime(primer_cand, native_3_prime_pos, side, seq_id, anchor_kmers, sequences_by_id, cfg):
    """Validates the primer against the k-mer map to ensure no high-Tm off-targets."""
    anchor_3 = primer_cand[-cfg['primer_3prime_anchor']:]
    
    for anchor_seq_id, anchor_pos, anchor_side in anchor_kmers.get(anchor_3, []):
        # Allow binding to the intended on-target site
        if anchor_side == side and anchor_seq_id == seq_id and anchor_pos == native_3_prime_pos:
            continue
            
        anchor_seq = str(sequences_by_id[anchor_seq_id].seq).upper()
        
        # Extract the off-target sequence depending on strand orientation
        if anchor_side == "forward":
            pos_misprime = anchor_seq[max(0, anchor_pos - cfg['max_primer_length'] + 1):anchor_pos + 1]
        else:
            pos_misprime = reverse_complement(anchor_seq[anchor_pos:min(len(anchor_seq), anchor_pos + cfg['max_primer_length'])])
            
        if primer3.calc_heterodimer_tm(primer_cand, reverse_complement(pos_misprime)) > cfg['max_misprime_tm']:
            return True
            
    return False

def check_local_misprime(primer_cand, seq_str, native_3_prime_pos, side, seq_id, anchor_kmers, cfg):
    """Checks for potential mispriming within the same potential fragment."""
    anchor_3 = primer_cand[-cfg['primer_3prime_anchor']:]
    
    for anchor_seq_id, anchor_pos, anchor_side in anchor_kmers.get(anchor_3, []):
        if anchor_seq_id != seq_id:
            continue
        if anchor_side != side:
            continue
        if anchor_pos == native_3_prime_pos:
            continue
        if anchor_pos - native_3_prime_pos > cfg['max_frag_length'] or native_3_prime_pos > anchor_pos:
            continue

        pos_misprime = seq_str[max(0, anchor_pos - cfg['max_primer_length'] + 1):anchor_pos + 1]
            
        if primer3.calc_heterodimer_tm(primer_cand, reverse_complement(pos_misprime)) > cfg['max_misprime_tm']:
            return True            
    return False

def find_primers_in_regions(seq_record, regions, side, anchor_kmers, sequences_by_id, cfg):
    """Finds primer candidates for a given sequence, regions, and orientation."""
    primer_candidates = []
    seq_id = seq_record.id
    original_seq = str(seq_record.seq).upper()
    
    if side == "forward":
        search_seq = original_seq
        search_regions = regions
    elif side == "reverse":
        search_seq = reverse_complement(original_seq)
        seq_len = len(original_seq)
        search_regions = [(seq_len - r[1], seq_len - r[0]) for r in regions]
    else:
        raise ValueError("Side must be 'forward' or 'reverse'")

    for region in search_regions:
        clamp_matches = cfg['clamp_pattern'].finditer(search_seq, region[0], region[1])
        
        for m in clamp_matches:
            range_start = max(m.start() - (cfg['max_primer_length'] - cfg['clamp_length']), region[0])
            range_end = m.start() - cfg['min_primer_length'] + cfg['clamp_length'] + 1
            
            # Extend 5' -> 3'
            for pr_start in range(range_end - 1, range_start - 1, -1):
                primer_cand = search_seq[pr_start:m.end()]
                if cfg['poly_x_pattern'].search(primer_cand): break
                if not cfg['five_prime_clamp_pattern'].search(primer_cand): continue
                
                gc_content = (primer_cand.count('G') + primer_cand.count('C')) / len(primer_cand)
                if gc_content < cfg['min_gc'] or gc_content > cfg['max_gc']: continue
                
                tm = primer3.calc_tm(primer_cand, **cfg['tm_params'])
                if tm > cfg['max_tm']: break
                if tm < cfg['min_tm']: continue
                if primer3.calc_hairpin_tm(primer_cand) > cfg['max_hairpin_tm']: continue
                if primer3.calc_homodimer_tm(primer_cand) > cfg['max_homodimer_tm']: continue
                if primer3.calc_end_stability(primer_cand, primer_cand).tm > cfg['max_3_self_tm']: continue
                
                if side == "forward":
                    native_3_prime_pos = m.end() - 1
                    pos_tuple = (pr_start, m.end())
                else:
                    native_3_prime_pos = len(original_seq) - m.end()
                    pos_tuple = (len(original_seq) - m.end(), len(original_seq) - pr_start)
                
                if check_local_misprime(primer_cand, search_seq, native_3_prime_pos, side, seq_id, anchor_kmers, cfg):
                    break
                    
                primer_candidates.append({
                    'seq': primer_cand,
                    'gc_content': gc_content,
                    'tm': tm,
                    'pos': pos_tuple,
                    'side': side
                })
                
    return primer_candidates

def find_primer_flanked_overlaps(seq_record, primer_regions, overlap_regions, anchor_kmers, sequences_by_id, cfg):
    """Finds all valid primer pairs that flank overlaps within the specified regions."""
    forward_primers = find_primers_in_regions(seq_record, primer_regions, "forward", anchor_kmers, sequences_by_id, cfg)
    reverse_primers = find_primers_in_regions(seq_record, primer_regions, "reverse", anchor_kmers, sequences_by_id, cfg)
    
    primer_flanked_overlaps = []

    for reg in overlap_regions:
        reg_forward_primers = list(filter(lambda x: x['pos'][0] >= reg[0] and x['pos'][1] <= reg[1], forward_primers))
        reg_reverse_primers = list(filter(lambda x: x['pos'][0] >= reg[0] and x['pos'][1] <= reg[1], reverse_primers))
        if len(reg_forward_primers) > 0 and len(reg_reverse_primers) > 0:
            for fwd in reg_forward_primers:
                for rev in reg_reverse_primers:
                    overlap_start = fwd['pos'][0]
                    overlap_end = rev['pos'][1]
                    if overlap_end - overlap_start >= cfg['min_overlap'] and overlap_end - overlap_start <= cfg['max_overlap']:
                        primer_flanked_overlaps.append({
                            'forward': fwd,
                            'reverse': rev,
                            'pos': (overlap_start, overlap_end)
                        })
    return primer_flanked_overlaps

# %%
anchor_kmers = build_anchor_kmers(sequences, CONFIG['primer_3prime_anchor'])

primer_flanked_overlaps = {}

for seq in sequences:
    print(seq.id)
    primer_flanked_overlaps[seq.id] = find_primer_flanked_overlaps(seq, bg_regions[seq.id], bg_regions[seq.id], anchor_kmers, sequences_by_id, CONFIG)

# %%
homfree_ranges = {}
for seq in sequences:
    seq_str = str(seq.seq).upper()
    homfree_ranges[seq.id] = compute_homfree_ranges(seq_str, CONFIG['kmer_size'])

# %%
def get_overlap_range(start, end, homfree_ranges, kmer_size):
    ov_end = min([el[1] for el in homfree_ranges[start:end - kmer_size + 1]])
    ov_start = max([el[0] for el in homfree_ranges[start + kmer_size - 1:end]])
    return ov_start, ov_end

for seq in sequences:
    seq_id = seq.id
    for ov in primer_flanked_overlaps[seq_id]:
        ov["range"] = get_overlap_range(ov["pos"][0], ov["pos"][1], homfree_ranges[seq_id], CONFIG['kmer_size'])
    primer_flanked_overlaps[seq_id] = [ov for ov in primer_flanked_overlaps[seq_id] if not( ov["range"][0] > ov["pos"][1] - CONFIG['min_frag_length'] or ov["range"][1] < ov["pos"][0] + CONFIG['min_frag_length'])]
# %%
# to bridge the complex repetitive gap we add two sort of primer-flanked, shorter overlaps
# only one of them is expected to be used, but let algorithm decide, which one
seq_id = 'Pereira_COS-SynDNA-f2'

pos = (16411, 16411 + 43)
seq_str = str(sequences_by_id[seq_id].seq)
print(seq_str[pos[0]:pos[1]])
primer_cand = "GTAGTCACATACCTGAAGAGGCAC"
print("forward one:")
print(primer3.calc_tm(primer_cand, **CONFIG['tm_params']))
print(primer3.calc_hairpin_tm(primer_cand))
print(primer3.calc_homodimer_tm(primer_cand))
print(primer3.calc_end_stability(primer_cand, primer_cand).tm)
# Well, that's a very nasty hairpin, to be honest...
forward_pr = {"seq": primer_cand,
                "gc_content": (primer_cand.count('G') + primer_cand.count('C')) / len(primer_cand),
                "tm": primer3.calc_tm(primer_cand, **CONFIG['tm_params']),
                "pos": (pos[0], pos[0] + len(primer_cand)),
                "side": "forward"}

primer_cand = reverse_complement("GGCAGAAAGTTTCACCTGTTCT")
print("reverse one:")
print(primer3.calc_tm(primer_cand, **CONFIG['tm_params']))
print(primer3.calc_hairpin_tm(primer_cand))
print(primer3.calc_homodimer_tm(primer_cand))
print(primer3.calc_end_stability(primer_cand, primer_cand).tm)
reverse_pr = {"seq": primer_cand,
                "gc_content": (primer_cand.count('G') + primer_cand.count('C')) / len(primer_cand),
                "tm": primer3.calc_tm(primer_cand, **CONFIG['tm_params']),
                "pos": (pos[1] - len(primer_cand), pos[1]),
                "side": "reverse"}
range_ov = get_overlap_range(pos[0], pos[1], homfree_ranges[seq_id], CONFIG['kmer_size'])
print(f"Overlap range: {range_ov}")
primer_flanked_overlaps[seq_id].append({
    "forward": forward_pr,
    "reverse": reverse_pr,
    "pos": (pos[0], pos[1]),
    "range": range_ov
})

pos = (16928, 16928 + 40)
seq_str = str(sequences_by_id[seq_id].seq)
print(seq_str[pos[0]:pos[1]])
primer_cand = "TTACTCACAACATACAGAGAAGCC"
print("forward two:")
print(primer3.calc_tm(primer_cand, **CONFIG['tm_params']))
print(primer3.calc_hairpin_tm(primer_cand))
print(primer3.calc_homodimer_tm(primer_cand))
print(primer3.calc_end_stability(primer_cand, primer_cand).tm)
# Well, that's a very nasty hairpin, to be honest...
forward_pr = {"seq": primer_cand,
                "gc_content": (primer_cand.count('G') + primer_cand.count('C')) / len(primer_cand),
                "tm": primer3.calc_tm(primer_cand, **CONFIG['tm_params']),
                "pos": (pos[0], pos[0] + len(primer_cand)),
                "side": "forward"}

primer_cand = reverse_complement("GAGAAGCCGAGTATTTTCTACCAA")
print("reverse two:")
print(primer3.calc_tm(primer_cand, **CONFIG['tm_params']))
print(primer3.calc_hairpin_tm(primer_cand))
print(primer3.calc_homodimer_tm(primer_cand))
print(primer3.calc_end_stability(primer_cand, primer_cand).tm)
reverse_pr = {"seq": primer_cand,
                "gc_content": (primer_cand.count('G') + primer_cand.count('C')) / len(primer_cand),
                "tm": primer3.calc_tm(primer_cand, **CONFIG['tm_params']),
                "pos": (pos[1] - len(primer_cand), pos[1]),
                "side": "reverse"}
range_ov = get_overlap_range(pos[0], pos[1], homfree_ranges[seq_id], CONFIG['kmer_size'])
print(f"Overlap range: {range_ov}")
primer_flanked_overlaps[seq_id].append({
    "forward": forward_pr,
    "reverse": reverse_pr,
    "pos": (pos[0], pos[1]),
    "range": range_ov
})
# %%
gap_threshold = 1.5 * CONFIG['max_frag_length']
gaps = {}

for seq in sequences:
    seq_len = len(seq.seq)
    gaps[seq.id] = []
    
    if primer_flanked_overlaps[seq.id]:
        overlaps = sorted([
            ov['pos'] 
            for ov in primer_flanked_overlaps[seq.id]
        ])
        
        if overlaps[0][0] > gap_threshold:
            gaps[seq.id].append((0, overlaps[0][0]))
        
        for i in range(len(overlaps) - 1):
            gap_start = overlaps[i][1]
            gap_end = overlaps[i + 1][0]
            gap_size = gap_end - gap_start
            if gap_size > gap_threshold:
                gaps[seq.id].append((gap_start, gap_end))
        
        if overlaps[-1][1] < seq_len - gap_threshold:
            gaps[seq.id].append((overlaps[-1][1], seq_len))
    else:
        gaps[seq.id].append((0, seq_len))

for seq in sequences:
    print(f"\n{seq.id}: {len(gaps[seq.id])} gaps found")
    for i, (gap_start, gap_end) in enumerate(gaps[seq.id]):
        print(f"  Gap {i+1}: {gap_start}-{gap_end} (length: {gap_end - gap_start})")

# %%
segment_overlaps = {}
for seq in sequences:
    print(seq.id)
    cur_overlaps = []
    for region in cur_seq_regions[seq.id]:
        if region[1] - region[0] < CONFIG['max_overlap'] and region[1] - region[0] >= CONFIG['min_overlap']:
            cur_overlaps.append(region)
            continue
        for start in range(region[0], region[1] - CONFIG['max_overlap'] + 1, CONFIG['min_step']):
            cur_overlaps.append((start, start + CONFIG['max_overlap']))
        for end in range(region[1], region[0] + CONFIG['max_overlap'] - 1, -CONFIG['min_step']):
            cur_overlaps.append((end - CONFIG['max_overlap'], end))
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

    if left_overlap_len == 0:
        left_overlap_len = opt_overlap_len
    if right_overlap_len == 0:
        right_overlap_len = opt_overlap_len

    penalty += abs(left_overlap_len - opt_overlap_len) / max_overlap_len_diff * 0.1
    penalty += abs(right_overlap_len - opt_overlap_len) / max_overlap_len_diff * 0.1

    return round(penalty * 100)

def get_edge_clamp(left_ov, right_ov):
    return left_ov['left_clamp'] + right_ov['right_clamp']


# %%
#offsets: for primers 40 in total, for segments 110 (from each side)

def get_all_overlaps(primer_flanked_overlaps, segment_overlaps, seq_len):
    all_overlaps = (
        [
            {'pos': ov, 'type': 'seg', 'range': (0, seq_len)} 
            for ov in segment_overlaps
        ] +
        [
            {
                'pos': ov['pos'], 
                'type': 'fr', 
                'pr_ind': i,
                'left_clamp': 3 if ov['forward']['seq'][0] == 'G' else 2,
                'right_clamp': 3 if ov['reverse']['seq'][0] == 'G' else 2,
                'pr_left_len': len(ov['forward']['seq']),
                'pr_right_len': len(ov['reverse']['seq']),
                'tm_left': ov['forward']['tm'],
                'tm_right': ov['reverse']['tm'],
                'range': ov['range']
            } 
            for i, ov in enumerate(primer_flanked_overlaps)
        ] +
        [
            {'pos': (0, 0), 'type': 'seg', 'range': (0, seq_len)}, 
            {'pos': (seq_len, seq_len), 'type': 'seg', 'range': (0, seq_len)}
        ]
    )
    return sorted(all_overlaps, key=lambda x: x['pos'][0])

def get_overlap_states(all_overlaps, min_length, max_length, min_seg_length, opt_seg_length, seq_len,
                       segment_offset_left, segment_offset_right, max_overlap, min_overlap,
                       seq_offset_left, seq_offset_right):

    # TO DO: move to config, but calculate somehow based on input pars
    opt_primer_len = 20
    max_primer_len_diff = 3

    overlap_states = [{} for _ in range(len(all_overlaps))]

    overlap_states[0][(0, seq_len)] = (0, -1, None) # state: (start_seg, end_range): (weight, predecessor_ind, predecessor_state)

    for i, ov in enumerate(all_overlaps):

        if i % 1000 == 0:
            print(f"Processing overlap {i}/{len(all_overlaps)}")
        if len(overlap_states[i]) == 0:
            continue
        j = i + 1
        too_far = False
        while j < len(all_overlaps) and not too_far:
            next_ov = all_overlaps[j]
            if next_ov['pos'][0] - max(ov['pos'][0], 0) > max_length:
                too_far = True
                continue

            length = next_ov['pos'][1] - max(ov['pos'][0], 0)
            if min_length <= length <= max_length:
                edge_weight = 100 + get_edge_penalty(ov, next_ov, 
                                                    opt_primer_len, max_primer_len_diff, 
                                                    max_overlap, max_overlap - min_overlap)
                # split in two: check if next_ov starts a new segment or if it continues the current one
                # segment start loop
                if next_ov['type'] == 'seg' or next_ov['range'] == (0, seq_len):
                    for state, (cur_weight, _, _) in overlap_states[i].items():
                        if next_ov['range'][0] > state[0] or next_ov['pos'][1] > state[1]:
                            continue
                        seg_len = next_ov['pos'][1] - state[0]
                        if seg_len < min_seg_length:
                            continue

                        #segments must have extra space for flanking stuff
                        extra_length = segment_offset_right
                        if state[0] == ov['pos'][0]:
                            extra_length += segment_offset_left
                        if ov['pos'] == (0, 0):
                            extra_length += seq_offset_left
                        if next_ov['pos'] == (seq_len, seq_len):
                            extra_length += seq_offset_right
                        if length + extra_length > max_length:
                            continue

                        extra_penalty = max(opt_seg_length - seg_len, 0) / opt_seg_length * 500

                        new_weight = cur_weight + edge_weight + extra_penalty
                        new_state = (next_ov['pos'][0], next_ov['range'][1])

                        existing = overlap_states[j].get(new_state)
                        if existing is None or new_weight < existing[0]:
                            overlap_states[j][new_state] = (new_weight, i, state)

                # continue segment loop
                if next_ov['type'] == 'fr':
                    for state, (cur_weight, _, _) in overlap_states[i].items():
                        if next_ov['range'][0] > state[0] or next_ov['pos'][1] > state[1]:
                            continue
                        if next_ov['pos'][1] - state[0] > opt_seg_length:
                            # in this case only start segment option is possible
                            continue
                        extra_length = 0
                        if ov['pos'] == (0, 0):
                            extra_length += seq_offset_left
                        if next_ov['pos'] == (seq_len, seq_len):
                            extra_length += seq_offset_right
                        if state[0] == ov['pos'][0]:
                            extra_length += segment_offset_left
                        
                        if length + extra_length > max_length:
                            # if we are at the start of the segment, we need extra space on the left
                            continue

                        new_weight = cur_weight + edge_weight
                        new_state = (state[0], min(next_ov['range'][1], state[1]))

                        existing = overlap_states[j].get(new_state)
                        if existing is None or new_weight < existing[0]:
                            overlap_states[j][new_state] = (new_weight, i, state)
                    
            j += 1
    return overlap_states

# %%
full_overlap_states = {}
full_all_overlaps  = {}
for seq in sequences:
    print(seq.id)
    full_all_overlaps[seq.id] = get_all_overlaps(primer_flanked_overlaps[seq.id], segment_overlaps[seq.id], len(str(seq.seq)))

    full_overlap_states[seq.id] = get_overlap_states(full_all_overlaps[seq.id], CONFIG['min_frag_length'], CONFIG['max_frag_length'], CONFIG['min_segment_length'], CONFIG['opt_segment_length'], len(str(seq.seq)),
                                                      CONFIG['segment_offset_left'], CONFIG['segment_offset_right'], CONFIG['max_overlap'], CONFIG['min_overlap'],
                                                      CONFIG['seq_offset_left'], CONFIG['seq_offset_right'])


# %%
for seq_id in full_overlap_states:
    print(seq_id)
    print(full_overlap_states[seq_id][-1])

# %%
def reconstruct_shortest_path(overlap_states, end_index=None):
    path = []
    current_index = len(overlap_states) - 1 if end_index is None else end_index
    current_state = min(overlap_states[current_index].items(), key=lambda x: x[1][0])[0]

    while current_index != -1:
        path.append((current_index, current_state))
        _, predecessor_index, predecessor_state = overlap_states[current_index][current_state]
        current_index = predecessor_index
        current_state = predecessor_state

    return path[::-1]

seq_id = 'Pereira_COS-SynDNA-f2'

#full_all_overlaps[seq_id] = get_all_overlaps(primer_flanked_overlaps[seq_id], segment_overlaps[seq_id], len(str(sequences_by_id[seq_id].seq)))
#
#full_overlap_states[seq_id] = get_overlap_states(full_all_overlaps[seq_id], CONFIG['min_frag_length'], CONFIG['max_frag_length'], CONFIG['min_segment_length'], CONFIG['opt_segment_length'], len(str(sequences_by_id[seq_id].seq)),
#                                                    CONFIG['segment_offset_left'], CONFIG['segment_offset_right'], CONFIG['max_overlap'], CONFIG['min_overlap'] )


last_index = 0
for i in range(len(full_overlap_states[seq_id])):
    if len(full_overlap_states[seq_id][i]) > 0:
        last_index = i
reconstruct_shortest_path(full_overlap_states[seq_id], end_index=last_index)

# %%
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os


def _is_virtual_sentinel(ov, seq_len):
    """True for the two dummy bookend overlaps added in get_all_overlaps."""
    return ov['pos'][1] <= 0 or ov['pos'][0] >= seq_len


def generate_annotated_genbank(seq_id, sequences_by_id, full_overlap_states, full_all_overlaps,
                                cur_seq_regions, primer_flanked_overlaps_dict, end_index=None):
    """Generate an annotated GenBank record for a sequence with a valid assembly path.

    Features annotated:
      - misc_feature : non-homologous regions (cur_seq_regions)
      - repeat_region: every real overlap node in the path (excl. virtual sentinels)
      - gene         : each assembly segment (between n_frags==0 nodes)
      - misc_feature : each PCR fragment within a segment
      - primer_bind  : fwd/rev primers for 'fr'-type overlaps strictly inside a segment
    """
    states = full_overlap_states[seq_id]
    all_overlaps = full_all_overlaps[seq_id]

    check_index = end_index if end_index is not None else len(states) - 1
    if not states[check_index]:
        return None

    path = reconstruct_shortest_path(states, end_index)
    seq_record = sequences_by_id[seq_id]
    seq_len = len(seq_record.seq)

    new_record = SeqRecord(
        seq_record.seq,
        id=seq_id[:16],
        name=seq_id[:16],
        description=seq_id,
        annotations={"molecule_type": "DNA"},
    )

    features = []

    # --- Non-homologous regions ---
    for start, end in cur_seq_regions.get(seq_id, []):
        features.append(SeqFeature(
            FeatureLocation(start, end),
            type="misc_feature",
            qualifiers={"label": ["non_homologous"],
                        "note": ["non-homologous region"]},
        ))

    # --- All real overlap nodes ---
    for ov_idx, state in path:
        ov = all_overlaps[ov_idx]
        if _is_virtual_sentinel(ov, seq_len):
            continue
        ov_type = ov['type']
        label_prefix = "junc" if state[0] == ov['pos'][0] else "ov"
        features.append(SeqFeature(
            FeatureLocation(ov['pos'][0], ov['pos'][1]),
            type="repeat_region",
            qualifiers={
                "label":    [f"{label_prefix}_{ov['pos'][0]}_{ov['pos'][1]}"],
                "note":     [f"overlap type={ov_type}, seg_start={state[0]}, clamp={state[1]}"],
                "rpt_type": ["overlap"],
            },
        ))

    # --- Segments: slices of path between segment-start nodes ---
    # A node is a segment start when state[0] (seg_start) equals the overlap's own position
    boundary_path_indices = [
        i for i, (ov_idx, state) in enumerate(path)
        if state[0] == all_overlaps[ov_idx]['pos'][0]
    ]
    # For partial paths that end mid-segment, include the last node as a boundary
    if path:
        last_ov_idx, last_state = path[-1]
        if last_state[0] != all_overlaps[last_ov_idx]['pos'][0]:
            boundary_path_indices.append(len(path) - 1)

    segment_num = 0
    for si in range(len(boundary_path_indices) - 1):
        pa = boundary_path_indices[si]
        pb = boundary_path_indices[si + 1]
        seg_slice = path[pa : pb + 1]
        # seg_slice[0]  = left junction (n_frags==0)
        # seg_slice[-1] = right junction (n_frags==0)
        # seg_slice[1:-1] = interior overlaps of this segment

        left_ov  = all_overlaps[seg_slice[0][0]]
        right_ov = all_overlaps[seg_slice[-1][0]]

        # Skip vacuous slices produced entirely from virtual sentinels
        if _is_virtual_sentinel(left_ov, seq_len) and _is_virtual_sentinel(right_ov, seq_len):
            continue

        seg_start = max(left_ov['pos'][0], 0)
        seg_end   = min(right_ov['pos'][1], seq_len)

        segment_num += 1
        n_frags = len(seg_slice) - 1

        # Segment feature
        features.append(SeqFeature(
            FeatureLocation(seg_start, seg_end),
            type="gene",
            qualifiers={
                "label": [f"Segment_{segment_num}"],
                "note":  [f"Assembly segment {segment_num}, {n_frags} fragment(s)"],
            },
        ))

        # Fragment features (one per edge in the slice)
        for j in range(1, len(seg_slice)):
            l_ov = all_overlaps[seg_slice[j - 1][0]]
            r_ov = all_overlaps[seg_slice[j][0]]
            frag_start = max(l_ov['pos'][0], 0)
            frag_end   = min(r_ov['pos'][1], seq_len)
            features.append(SeqFeature(
                FeatureLocation(frag_start, frag_end),
                type="misc_feature",
                qualifiers={
                    "label": [f"S{segment_num}_F{j}"],
                    "note":  [f"Segment {segment_num}, fragment {j}"],
                },
            ))

        # Primer features for interior 'fr'-type overlaps
        for j in range(1, len(seg_slice) - 1):
            ov_idx, _ = seg_slice[j]
            ov = all_overlaps[ov_idx]
            if ov['type'] != 'fr':
                continue

            pr_data = primer_flanked_overlaps_dict[seq_id][ov['pr_ind']]
            fwd = pr_data['forward']
            rev = pr_data['reverse']

            features.append(SeqFeature(
                FeatureLocation(fwd['pos'][0], fwd['pos'][1], strand=1),
                type="primer_bind",
                qualifiers={
                    "label":    [f"S{segment_num}_Ov{j}_F"],
                    "note":     [f"Fwd primer Tm={fwd['tm']:.1f} GC={fwd['gc_content']:.2f}"],
                    "sequence": [fwd['seq']],
                },
            ))
            features.append(SeqFeature(
                FeatureLocation(rev['pos'][0], rev['pos'][1], strand=-1),
                type="primer_bind",
                qualifiers={
                    "label":    [f"S{segment_num}_Ov{j}_R"],
                    "note":     [f"Rev primer Tm={rev['tm']:.1f} GC={rev['gc_content']:.2f}"],
                    "sequence": [rev['seq']],
                },
            ))

    new_record.features = sorted(features, key=lambda f: f.location.start)
    return new_record


def extract_fragment_records(seq_id, sequences_by_id, path, all_overlaps):
    """Return a list of SeqRecords, one per fragment edge in the path."""
    seq_record = sequences_by_id[seq_id]
    seq_len = len(seq_record.seq)
    frag_records = []
    frag_num = 0
    for i in range(1, len(path)):
        prev_ov = all_overlaps[path[i - 1][0]]
        curr_ov = all_overlaps[path[i][0]]
        frag_start = max(prev_ov['pos'][0], 0)
        frag_end   = min(curr_ov['pos'][1], seq_len)
        if frag_start >= frag_end:
            continue
        frag_num += 1
        frag_records.append(SeqRecord(
            seq_record.seq[frag_start:frag_end],
            id=f"{seq_id}_F{frag_num}",
            name=f"{seq_id}_F{frag_num}"[:16],
            description=f"fragment {frag_num} pos={frag_start}-{frag_end}",
        ))
    return frag_records


# Write GenBank and FASTA files for every sequence
os.makedirs("annotated", exist_ok=True)

for seq_id in full_overlap_states:
    states = full_overlap_states[seq_id]
    end_index = None
    partial = False

    if not states[-1]:
        # Find the last reachable overlap index
        last_idx = None
        for idx in range(len(states) - 1, -1, -1):
            if states[idx]:
                last_idx = idx
                break
        if last_idx is None or last_idx == 0:
            print(f"{seq_id}: no solution found, skipping")
            continue
        end_index = last_idx
        partial = True
        print(f"{seq_id}: partial solution found (up to overlap index {last_idx})")

    record = generate_annotated_genbank(
        seq_id,
        sequences_by_id,
        full_overlap_states,
        full_all_overlaps,
        cur_seq_regions,
        primer_flanked_overlaps,
        end_index=end_index,
    )
    if record is None:
        print(f"{seq_id}: no solution found, skipping")
        continue

    suffix = "_partial" if partial else ""

    # Write GenBank
    out_gb = f"annotated/{seq_id}{suffix}.gb"
    with open(out_gb, "w") as fh:
        SeqIO.write(record, fh, "genbank")
    print(f"Written: {out_gb}  ({len(record.features)} features)")

    # Write FASTA fragments
    path = reconstruct_shortest_path(states, end_index)
    frag_records = extract_fragment_records(seq_id, sequences_by_id, path, full_all_overlaps[seq_id])
    out_fa = f"annotated/{seq_id}{suffix}_fragments.fa"
    with open(out_fa, "w") as fh:
        SeqIO.write(frag_records, fh, "fasta")
    print(f"Written: {out_fa}  ({len(frag_records)} fragments)")


# %%
list(filter(lambda ov: ov['pos'][0] > 15000 and ov['pos'][1] < 17000 and ov['range'][1] > 17000, primer_flanked_overlaps[seq_id]))

seq_id = "Pereira_COS-SynDNA-f2"
bridge = [False] * len(str(sequences_by_id[seq_id].seq))
for i, range in enumerate(homfree_ranges[seq_id]):
    bridge[i] = range[0] < 15500 and range[1] > 17300


bridge_slice = bridge[16000:17000]
consecutive_true = []
start = None
for i, val in enumerate(bridge_slice):
    if val:
        if start is None:
            start = i
    else:
        if start is not None:
            consecutive_true.append({'position': 16000 + start, 'size': i - start})
            start = None
if start is not None:
    consecutive_true.append({'position': 16000 + start, 'size': len(bridge_slice) - start})

for stretch in consecutive_true:
    print(f"Position: {stretch['position']}, Size: {stretch['size']}")

# %%
primer_cand = "TTACTCACAACATACAGAGAAGCC"
print(primer3.calc_tm(primer_cand, **CONFIG['tm_params']))
print(primer3.calc_hairpin_tm(primer_cand))
print(primer3.calc_homodimer_tm(primer_cand))
print(primer3.calc_end_stability(primer_cand, primer_cand).tm)