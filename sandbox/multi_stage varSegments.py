# %%
from Bio import SeqIO
import primer3
from choppy.homology_finder import load_trie, create_kmer_trie, merge_tries, find_non_homologous_regions, find_local_non_homologous_regions, get_region_intersects
from choppy.primer_flanked_overlaps import get_primer_overlaps_from_seq, check_primer
import re
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib.pyplot as plt
# %%
all_sequences = list(SeqIO.parse("data/20240414-forSveta.fa", "fasta"))

seq_names = ["Maizel_COS-AT1G27430-GYF2", "Maizel_COS-SETH5", 
             "Pereira_COS-SynDNA-f1", "Pereira_COS-SynDNA-f2",
             "PV252688r_p6utr", "PQ537341r_p6utr,"
             "PQ488556r_p6utr", "PX021458r_p6utr",
             "MZ289137_rep", "OR500095r_naive"]

sequences = [seq for seq in all_sequences if seq.id in seq_names]

# %%
min_step = 10
min_overlap = 50
max_overlap = 100

min_length = 100 - 40
min_length = max(min_length, min_overlap*2 + 40)
max_length = 1000 - 40
n_frags_max = 5
n_frags_min = 2
kmer_size = 15
primer_3prime_anchor = 6

threshold = 50
pr_threshold = 17
neighbourhood_size = 4500
# %%
for seq in sequences:
    seq.seq = seq.seq.upper()

bg_trie = load_trie("../data/S_cerevisiae-R64-GCA_000146045_cat_15.marisa")
full_seq_trie = create_kmer_trie(sequences, kmer_size, bg=False)

seq_tries = {}
for seq in sequences:
    print(f"Processing {seq.id}...")
    trie = create_kmer_trie(seq, kmer_size, bg=False)
    seq_tries[seq.id] = trie 
# %%
bg_regions = {}
full_seq_regions = {}
cur_seq_regions = {}
local_neigh_regions = {}
for seq in sequences:
    bg_regions[seq.id] = find_non_homologous_regions(seq, bg_trie, [], kmer_size, threshold=threshold)
    full_seq_regions[seq.id] = find_non_homologous_regions(seq, full_seq_trie, bg_trie, kmer_size, threshold=pr_threshold)
    cur_seq_regions[seq.id] = find_non_homologous_regions(seq, seq_tries[seq.id], bg_trie, kmer_size, threshold=threshold)
    local_neigh_regions[seq.id] = find_local_non_homologous_regions(seq, kmer_size, threshold, neighbourhood_size)
    local_neigh_regions[seq.id] = get_region_intersects(local_neigh_regions[seq.id], bg_regions[seq.id], threshold)

# %%
CONFIG = {
    'min_primer_length': 17,
    'max_primer_length': 23,
    'primer_3prime_anchor': 6,
    'clamp_length': 3,
    'min_gc': 0.3,
    'max_gc': 0.7,
    'min_tm': 57.0,
    'max_tm': 62.0,
    'max_hairpin_tm': 24.0,
    'max_homodimer_tm': 45.0,
    'max_3_self_tm': 35.0,
    'max_misprime_tm': 47.0,
    'min_overlap': 50,
    'max_overlap': 100,
    'poly_x_pattern': re.compile(r'(A{5,}|T{5,}|G{5,}|C{5,})'),
    'clamp_pattern': re.compile(r'[GC][AT][GC]|[AT][GC][GC]'),
    'five_prime_clamp_pattern': re.compile(r'^[GC]')
}
# %%
def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans('ATGC', 'TACG'))

def build_anchor_kmers(sequences, anchor_len):
    """Builds the 3' k-mer hash map for fast off-target screening."""
    anchor_kmers = defaultdict(list)
    for seq in sequences:
        seq_str = str(seq.seq).upper()
        for pos in range(len(seq_str) - anchor_len + 1):
            kmer_fwd = seq_str[pos:pos + anchor_len]
            kmer_rev = reverse_complement(kmer_fwd)
            
            anchor_kmers[kmer_fwd].append((seq.id, pos + anchor_len - 1, "forward"))
            anchor_kmers[kmer_rev].append((seq.id, pos, "reverse"))
            
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
            
        if primer3.calc_heterodimer_tm(primer_cand, pos_misprime) > cfg['max_misprime_tm']:
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
                
                tm = primer3.calc_tm(primer_cand)
                if tm > cfg['max_tm']: break
                if tm < cfg['min_tm']: continue
                if primer3.calc_hairpin_tm(primer_cand) > cfg['max_hairpin_tm']: break
                if primer3.calc_homodimer_tm(primer_cand) > cfg['max_homodimer_tm']: break
                if primer3.calc_end_stability(primer_cand, primer_cand).tm > cfg['max_3_self_tm']: break
                
                if side == "forward":
                    native_3_prime_pos = m.end() - 1
                    pos_tuple = (pr_start, m.end())
                else:
                    native_3_prime_pos = len(original_seq) - m.end()
                    pos_tuple = (len(original_seq) - m.end(), len(original_seq) - pr_start)
                
                if check_misprime(primer_cand, native_3_prime_pos, side, seq_id, anchor_kmers, sequences_by_id, cfg):
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

# %% Reverse check: brute-force misprime Tm verification
sequences_by_id = {seq.id: seq for seq in sequences}
anchor_kmers = build_anchor_kmers(sequences, CONFIG['primer_3prime_anchor'])
# %%
primer_flanked_overlaps = {}

for seq in sequences:
    print(seq.id)
    primer_flanked_overlaps[seq.id] = find_primer_flanked_overlaps(seq, local_neigh_regions[seq.id], full_seq_regions[seq.id], anchor_kmers, sequences_by_id, CONFIG)

# %%
for ov in primer_flanked_overlaps[sequences[0].id]:
    print(ov["reverse"]["seq"])

# %% Identify gaps in primer_flanked_overlaps coverage
gap_threshold = 1.5 * max_length
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
RELAXED_CONFIG = {
    'min_primer_length': 17,
    'max_primer_length': 28,
    'primer_3prime_anchor': 6,
    'clamp_length': 1,
    'min_gc': 0.3,
    'max_gc': 0.7,
    'min_tm': 57.0,
    'max_tm': 65.0,
    'max_hairpin_tm': 24.0,
    'max_homodimer_tm': 45.0,
    'max_3_self_tm': 35.0,
    'max_misprime_tm': 50.0,
    'min_overlap': 50,
    'max_overlap': 120,
    'poly_x_pattern': re.compile(r'(A{5,}|T{5,}|G{5,}|C{5,})'),
    'clamp_pattern': re.compile(r'[GC]'),
    'five_prime_clamp_pattern': re.compile(r'^[ATGC]')
}

primer_flanked_overlaps_relaxed = {}

for seq in sequences:
    print(seq.id)
    if len(gaps[seq.id]) > 0:
        primer_regions = list(filter(lambda x: any(gap_start <= x[0] < gap_end or gap_start < x[1] <= gap_end for gap_start, gap_end in gaps[seq.id]), local_neigh_regions[seq.id]))
        overlap_regions = list(filter(lambda x: any(gap_start <= x[0] < gap_end or gap_start < x[1] <= gap_end for gap_start, gap_end in gaps[seq.id]), local_neigh_regions[seq.id]))
        primer_flanked_overlaps_relaxed[seq.id] = find_primer_flanked_overlaps(seq, primer_regions, overlap_regions, anchor_kmers, sequences_by_id, RELAXED_CONFIG)
        primer_flanked_overlaps_relaxed[seq.id] = [ov for ov in primer_flanked_overlaps_relaxed[seq.id] 
                                                   if any(gap_start <= ov['pos'][0] < gap_end and gap_start < ov['pos'][1] <= gap_end 
                                                          for gap_start, gap_end in gaps[seq.id])]
# %%
gap_threshold = 1.5 * max_length
gaps = {}

for seq in sequences:
    seq_len = len(seq.seq)
    gaps[seq.id] = []
    
    if primer_flanked_overlaps[seq.id]:
        overlaps = sorted([
            ov['pos'] 
            for ov in primer_flanked_overlaps[seq.id] + primer_flanked_overlaps_relaxed.get(seq.id, [])
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
#gc_clamp_pattern_relaxed = re.compile(r'[GC]')
gc_clamp_pattern_relaxed = re.compile(r'[GCAT]')

primer_flanked_overlaps_relaxed = {}
for seq in sequences:
    if len(gaps[seq.id]) > 0:
        candidate_regions = list(filter(lambda x: any(gap_start <= x[0] < gap_end or gap_start < x[1] <= gap_end for gap_start, gap_end in gaps[seq.id]), bg_regions[seq.id]))
        print(candidate_regions)        
        overlaps_relaxed = get_primer_overlaps_from_seq(str(seq.seq), candidate_regions, 
                                                        gc_clamp_pattern_relaxed, gc_clamp_pattern_relaxed, 1,
                                                        five_prime_clamp_pattern_forward=None, max_heterodimer_tm=90.0, 
                                                        five_prime_clamp_pattern_reverse=None, min_overlap=35, max_homodimer_tm=60.0,
                                                        max_tm=65.0, 
                                                        min_length=min_primer_length, max_length=max_primer_length+10)
        
        
        primer_flanked_overlaps_relaxed[seq.id] = [ov for ov in overlaps_relaxed 
                           if any(gap_start <= ov['overlap_start'] < gap_end and gap_start < ov['overlap_end'] <= gap_end 
                                  for gap_start, gap_end in gaps[seq.id])]       

# %%

# %% Dot plot: Pereira_COS-SynDNA-f2 gap region vs full sequence (10-mer matches)
_dotplot_seq_id = "PX021458r_p6utr"
_seq = next(s for s in sequences if s.id == _dotplot_seq_id)
_full_seq = str(_seq.seq)
_gap_start, _gap_end = gaps[_dotplot_seq_id][0]
_gap_seq = _full_seq[_gap_start:_gap_end]

_k = 15

# Index all k-mers in the gap sequence
_gap_kmer_positions = {}
for _j in range(len(_gap_seq) - _k + 1):
    _km = _gap_seq[_j:_j + _k]
    _gap_kmer_positions.setdefault(_km, []).append(_j)

_comp_table = str.maketrans('ACGT', 'TGCA')
def _revcomp(s):
    return s.translate(_comp_table)[::-1]

_fw_x, _fw_y = [], []
_rc_x, _rc_y = [], []

for _i in range(len(_full_seq) - _k + 1):
    _km = _full_seq[_i:_i + _k]
    if _km in _gap_kmer_positions:
        for _j in _gap_kmer_positions[_km]:
            _fw_x.append(_i)
            _fw_y.append(_j)
    _km_rc = _revcomp(_km)
    if _km_rc in _gap_kmer_positions:
        for _j in _gap_kmer_positions[_km_rc]:
            _rc_x.append(_i)
            _rc_y.append(_j)

fig, ax = plt.subplots(figsize=(14, 5))
ax.scatter(_fw_x, _fw_y, s=1, c='steelblue', label='forward', alpha=0.6, linewidths=0)
ax.scatter(_rc_x, _rc_y, s=1, c='tomato', label='revcomp', alpha=0.6, linewidths=0)

# Overlay regions from cur_seq_regions and bg_regions that intersect the gap (y-axis = gap-relative coords)
_cur_in_gap = [(max(r[0], _gap_start) - _gap_start, min(r[1], _gap_end) - _gap_start)
               for r in local_neigh_regions[_dotplot_seq_id]
               if r[0] < _gap_end and r[1] > _gap_start]
_bg_in_gap  = [(max(r[0], _gap_start) - _gap_start, min(r[1], _gap_end) - _gap_start)
               for r in bg_regions[_dotplot_seq_id]
               if r[0] < _gap_end and r[1] > _gap_start]

_cur_label_done = False
_bg_label_done  = False
for _ys, _ye in _cur_in_gap:
    ax.axhspan(_ys, _ye, color='limegreen', alpha=0.20,
               label='cur_seq_regions' if not _cur_label_done else None)
    _cur_label_done = True
for _ys, _ye in _bg_in_gap:
    ax.axhspan(_ys, _ye, color='orange', alpha=0.20,
               label='bg_regions' if not _bg_label_done else None)
    _bg_label_done = True

ax.set_xlabel(f"Position in {_dotplot_seq_id} (full sequence, bp)")
ax.set_ylabel(f"Position in gap region\n(gap {_gap_start}–{_gap_end}, relative, bp)")
ax.set_title(f"Dot plot: gap region vs full sequence — {_k}-mer matches")
ax.legend(markerscale=6)
plt.tight_layout()
plt.show()
# %%
cur_region = 0
for i, x in enumerate(_fw_x):
    if x == _fw_y[i] + _gap_start:
        continue
    while cur_seq_regions[_dotplot_seq_id][cur_region][1] < x:
        cur_region += 1
    if cur_seq_regions[_dotplot_seq_id][cur_region][0] <= x < cur_seq_regions[_dotplot_seq_id][cur_region][1] - 14:
        print(f"Forward match at {_fw_x[i]}/{_fw_y[i] + _gap_start} is in cur_seq_region {cur_seq_regions[_dotplot_seq_id][cur_region]}")



# %%
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


# %%

segment_offset_left = 78
segment_offset_right = 108
#offsets: for primers 40 in total, for segments 110 (from each side)

def get_all_overlaps(primer_flanked_overlaps, segment_overlaps, seq_len):
    all_overlaps = (
        [
            {'pos': ov, 'type': 'seg'} 
            for ov in segment_overlaps
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
            for i, ov in enumerate(primer_flanked_overlaps)
        ] +
        [
            {'pos': (-100, 0), 'type': 'seg'}, 
            {'pos': (seq_len, seq_len + 100), 'type': 'seg'}
        ]
    )
    return sorted(all_overlaps, key=lambda x: x['pos'][0])

def get_overlap_states(all_overlaps, min_length, max_length, n_frags_min, n_frags_max,
                       segment_offset_left, segment_offset_right, max_overlap, min_overlap, 
                       max_primer_length, min_primer_length):

    opt_primer_len = (max_primer_length + min_primer_length) / 2
    max_primer_len_diff = max_primer_length - min_primer_length

    overlap_states = [{} for _ in range(len(all_overlaps))]

    overlap_states[0][(0, 0)] = (0, -1, None) # state: (n_frags, clamp_code): (weight, predecessor_ind, predecessor_state)

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
                for state, (cur_weight, _, _) in overlap_states[i].items():
                    new_segment_started = False
                    new_n_frag = state[0] + 1
                    if new_n_frag > n_frags_max:
                        new_n_frag = 0
                        new_segment_started = True

                    new_clamp = state[1]
                    new_weight = cur_weight + edge_weight

                    if next_ov['type'] == 'seg':
                        if new_n_frag >= n_frags_min:
                            new_n_frag = 0
                            new_segment_started = True
                        else:
                            continue
                    # clamp-based restrictions
                    elif new_n_frag == 1:
                        # we have just started
                        new_clamp = next_ov['right_clamp']
                    else:
                        if not new_segment_started:
                            # we first attempt to continue segment, and if not possible, stop it prematurely
                            failed = False
                            # in the middle of the segment, both clamps matter
                            edge_clamp = ov['left_clamp'] + next_ov['right_clamp']

                            # at this point we know that we are at a middle fragment, edge_clamp > 3
                            # if not, something went very wrong
                            if edge_clamp <= 3:
                                raise ValueError(f"Unexpected clamp code {edge_clamp} for edge between {ov} and {next_ov}")

                            if state[1] > 3 and state[1] != edge_clamp:
                                failed = True
                            else:
                                if state[1] == 2 and edge_clamp == 6: # C vs GG
                                    failed = True    
                                elif state[1] == 3 and edge_clamp == 4: # G vs CC
                                    failed = True
                            if failed and new_n_frag >= n_frags_min:
                                new_n_frag = 0
                                new_segment_started = True
                                new_weight += (n_frags_max - new_n_frag) * 100 # penalty for stopping segment prematurely
                
                            if new_clamp <= 3 and edge_clamp > 3:
                                new_clamp = edge_clamp

                            #technically, this should never happen, since it requires a very small segment
                            if new_clamp <= 3 and edge_clamp <= 3 and new_clamp != edge_clamp:
                                new_clamp += edge_clamp

                    if new_segment_started:
                        # end of the segmented = True, only left clamp of the current vertex matters (left of the new edge)
                        if ov['type'] == 'fr':
                            if ov['left_clamp'] == 2 and new_clamp == 6: # C vs GG
                                continue
                            elif ov['left_clamp'] == 3 and new_clamp == 4: # G vs CC
                                continue
                        new_clamp = 0
                    
                    # the first and last fragments in a segment require an extra offsest, therefore check once more for the length limit
                    if new_segment_started:
                        if length + segment_offset_left > max_length:
                            continue
                    if state[1] == 0:
                        if length + segment_offset_right > max_length:
                            continue
                    
                    new_state = (new_n_frag, new_clamp)
                    existing = overlap_states[j].get(new_state)
                    if existing is None or new_weight < existing[0]:
                        overlap_states[j][new_state] = (new_weight, i, state)
            j += 1
    return overlap_states


full_overlap_states = {}
full_all_overlaps  = {}
for seq in sequences:
    print(seq.id)
    full_all_overlaps[seq.id] = get_all_overlaps(primer_flanked_overlaps[seq.id], segment_overlaps[seq.id], len(str(seq.seq)))
    full_overlap_states[seq.id] = get_overlap_states(full_all_overlaps[seq.id], min_length, max_length, n_frags_min, n_frags_max,
                                                      segment_offset_left, segment_offset_right, max_overlap, min_overlap, 
                                                      max_primer_length, min_primer_length)

# %%
for seq_id in full_overlap_states:
    print(seq_id)
    print(full_overlap_states[seq_id][-1])

# %%
def reconstruct_shortest_path(overlap_states):
    path = []
    current_index = len(overlap_states) - 1
    current_state = min(overlap_states[current_index].items(), key=lambda x: x[1][0])[0]

    while current_index != -1:
        path.append((current_index, current_state))
        _, predecessor_index, predecessor_state = overlap_states[current_index][current_state]
        current_index = predecessor_index
        current_state = predecessor_state

    return path[::-1]

seq_id = "MZ289137_rep"
reconstruct_shortest_path(full_overlap_states[seq_id])

# %%
seq_id = "Pereira_COS-SynDNA-f2"
full_overlap_states[seq_id][-10:]

start = 15429
end = 17328

list(filter(lambda x: x[0] >= start and x[1] <= end, local_neigh_regions[seq_id]))
list(filter(lambda x: x[0] >= start and x[1] <= end, cur_seq_regions[seq_id]))
regs = list(filter(lambda x: x[0] >= start and x[1] <= end, primer_candidate_regions[seq_id]))

for reg in regs:
    print(sequences[3].seq[reg[0]:reg[1]])

# %%

last_reached_index = len(full_overlap_states[seq_id]) - 1
while last_reached_index >= 0 and len(full_overlap_states[seq_id][last_reached_index]) == 0:
    last_reached_index -= 1


# %%
total_frags = 0
for seq_id in full_overlap_states:
    if len(full_overlap_states[seq_id][-1]) > 0:
        path = reconstruct_shortest_path(full_overlap_states[seq_id])
        total_frags += len(path) - 1
print(total_frags)

# %%
sum([len(seq.seq) for seq in sequences])/ total_frags