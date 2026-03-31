import re
import primer3

calc_gc_content = lambda s: (s.count('G') + s.count('C')) / len(s)
get_poly_x_pattern = lambda max_poly_x_length: re.compile(
    r'(A{' + str(max_poly_x_length + 1) + r',}'
    r'|T{' + str(max_poly_x_length + 1) + r',}'
    r'|G{' + str(max_poly_x_length + 1) + r',}'
    r'|C{' + str(max_poly_x_length + 1) + r',})'
)
def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans('ATGC', 'TACG'))

def check_primer(
    primer, poly_x_pattern=None, poly_x_max_length=4,
    min_gc=0.3, max_gc=0.7, 
    min_tm=57, max_tm=62, 
    max_hairpin_tm=24, max_homodimer_tm=45, 
    max_3_self_tm=35.0, five_prime_clamp_pattern=None
):
    """ 
    Checks if a primer sequence meets certain criteria (checked in the order listed):
    - No homopolymer runs longer than poly_x_max_length
    - GC content between min_gc and max_gc
    - Melting temperature (Tm) between min_tm and max_tm
    - Hairpin Tm below max_hairpin_tm
    - Homodimer Tm below max_homodimer_tm
    - 3' self-complementarity Tm below max_3_self_tm
    - (optional) Presence of a 5' clamp sequence matching five_prime_clamp_pattern
    
    Returns a tuple (is_valid, gc_content, tm) where:
    - is_valid: True if the primer meets all criteria, False otherwise
    - gc_content: The calculated GC content of the primer (or None if invalid)
    - tm: The calculated melting temperature of the primer (or None if invalid)
    """

    if poly_x_pattern is None:
        poly_x_pattern = get_poly_x_pattern(poly_x_max_length)
        
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
    if five_prime_clamp_pattern is not None and not re.search(five_prime_clamp_pattern, primer):
        return False, gc_content, mt
    return True, gc_content, mt

def find_forward_primers(seq, allowed_regions, clamp_pattern, clamp_length,
                         max_length = 23, min_length = 18, max_poly_x_length=4,
                         min_gc=0.3, max_gc=0.7, min_tm=57, max_tm=62,
                         max_hairpin_tm=24, max_homodimer_tm=45, max_3_self_tm=35.0, 
                         five_prime_clamp_pattern=None):
    """
    Finds candidate forward primers in the given sequence that meet specified criteria.
    - seq: The input DNA sequence to search for primers.
    - allowed_regions: A list of (start, end) tuples specifying regions of the
      sequence to search for primers.
    - clamp_pattern: A regex pattern that identifies the clamp sequence in the
      primers.
    - clamp_length: The length of the clamp sequence that must be present in
      the primer.
    - max_length: The maximum allowed length of the primer.
    - min_length: The minimum allowed length of the primer.
    - max_poly_x_length: The maximum allowed length of homopolymer runs in the
      primer.
    - min_gc: The minimum allowed GC content of the primer.
    - max_gc: The maximum allowed GC content of the primer.
    - min_tm: The minimum allowed melting temperature (Tm) of the primer.
    - max_tm: The maximum allowed melting temperature (Tm) of the primer.
    - max_hairpin_tm: The maximum allowed hairpin Tm of the primer.
    - max_homodimer_tm: The maximum allowed homodimer Tm of the primer.
    - max_3_self_tm: The maximum allowed 3' self-complementarity Tm of the
      primer.

    Returns a list of candidate primers that meet the specified criteria, where
    each candidate is a dictionary with keys 'seq', 'gc_content', 'tm', and
    'pos' (the start and end positions of the primer in the input sequence).
    """
    forward_candidates = []
    poly_x_pattern = get_poly_x_pattern(max_poly_x_length)

    for region in allowed_regions:
        clamp_matches = clamp_pattern.finditer(seq, region[0], region[1])

        for m in clamp_matches:
            range_start = max(m.start() - (max_length - clamp_length), region[0])
            range_end = m.start() - min_length + clamp_length + 1
            for pr_start in range(range_start, range_end):
                primer_cand = seq[pr_start:m.end()]
                valid, gc_content, mt = check_primer(
                    primer_cand, poly_x_pattern,
                    poly_x_max_length=max_poly_x_length,
                    min_gc=min_gc, max_gc=max_gc,
                    min_tm=min_tm, max_tm=max_tm,
                    max_hairpin_tm=max_hairpin_tm,
                    max_homodimer_tm=max_homodimer_tm,
                    max_3_self_tm=max_3_self_tm,
                    five_prime_clamp_pattern=five_prime_clamp_pattern)
                if not valid:
                    continue
                forward_candidates.append({
                    'seq': primer_cand, 'gc_content': gc_content,
                    'tm': mt, 'pos': (pr_start, m.end())
                })

    return sorted(forward_candidates, key=lambda el: el["pos"][0])

def find_reverse_primers(seq, allowed_regions, clamp_pattern, clamp_length,
                         max_length = 23, min_length = 18, max_poly_x_length=4,
                         min_gc=0.3, max_gc=0.7, min_tm=57, max_tm=62,
                         max_hairpin_tm=24, max_homodimer_tm=45, max_3_self_tm=35.0,
                         five_prime_clamp_pattern=None):
    """
    Finds candidate reverse primers in the given sequence that meet specified criteria.
    - seq: The input DNA sequence to search for primers.
    - allowed_regions: A list of (start, end) tuples specifying regions of the
      sequence to search for primers.
    - clamp_pattern: A regex pattern that identifies the clamp sequence in the
      primers.
    - clamp_length: The length of the clamp sequence that must be present in
      the primer.
    - max_length: The maximum allowed length of the primer.
    - min_length: The minimum allowed length of the primer.
    - max_poly_x_length: The maximum allowed length of homopolymer runs in the
      primer.
    - min_gc: The minimum allowed GC content of the primer.
    - max_gc: The maximum allowed GC content of the primer.
    - min_tm: The minimum allowed melting temperature (Tm) of the primer.
    - max_tm: The maximum allowed melting temperature (Tm) of the primer.
    - max_hairpin_tm: The maximum allowed hairpin Tm of the primer.
    - max_homodimer_tm: The maximum allowed homodimer Tm of the primer.
    - max_3_self_tm: The maximum allowed 3' self-complementarity Tm of the
      primer.
    - five_prime_clamp_pattern: An optional regex pattern that identifies the
      5' clamp sequence in the primers.

    Returns a list of candidate primers that meet the specified criteria, where
    each candidate is a dictionary with keys 'seq', 'gc_content', 'tm', and
    'pos' (the start and end positions of the primer in the input sequence).
    """
    reverse_candidates = []
    poly_x_pattern = get_poly_x_pattern(max_poly_x_length)

    for region in allowed_regions:
        clamp_matches = clamp_pattern.finditer(seq, region[0], region[1])

        for m in clamp_matches:
            range_start = min(m.end() + (max_length - clamp_length), region[1])
            range_end = m.end() + min_length - clamp_length - 1
            for pr_end in range(range_start, range_end, -1):
                primer_cand = seq[m.start():pr_end]
                valid, gc_content, mt = check_primer(
                    primer_cand, poly_x_pattern,
                    poly_x_max_length=max_poly_x_length,
                    min_gc=min_gc, max_gc=max_gc,
                    min_tm=min_tm, max_tm=max_tm,
                    max_hairpin_tm=max_hairpin_tm,
                    max_homodimer_tm=max_homodimer_tm,
                    max_3_self_tm=max_3_self_tm,
                    five_prime_clamp_pattern=five_prime_clamp_pattern)
                if not valid:
                    continue
                reverse_candidates.append({
                    'seq': primer_cand, 'gc_content': gc_content,
                    'tm': mt, 'pos': (m.start(), pr_end)
                }) 

    return sorted(reverse_candidates, key=lambda el: el["pos"][1])

def get_primer_overlaps(forward_primers, reverse_primers, min_overlap=60, max_overlap=100, 
                        max_heterodimer_tm=45.0):
    """
    Identifies pairs of forward and reverse primers that have an overlap within a specified 
    range and do not form strong heterodimers.
    - forward_primers: A list of candidate forward primers, where each primer is a dictionary with keys 'seq', 'gc_content', 'tm', and 'pos' (the start and end positions of the primer in the input sequence).
    - reverse_primers: A list of candidate reverse primers, where each primer is a dictionary with keys 'seq', 'gc_content', 'tm', and 'pos' (the start and end positions of the primer in the input sequence).
    - min_overlap: The minimum allowed length of the overlap between the forward and reverse primers.
    - max_overlap: The maximum allowed length of the overlap between the forward and reverse primers.
    - max_heterodimer_tm: The maximum allowed melting temperature (Tm) of the heterodimer formed by the forward and reverse primers.    
    
    Returns a list of possible overlaps, where each overlap is a dictionary with keys 'forward', 'reverse', 'overlap_start', 'overlap_end', and 'overlap_length'. The 'forward' and 'reverse' keys contain the corresponding primer dictionaries from the input lists.  
    """

    possible_overlaps = []
    
    for i, primer1 in enumerate(forward_primers):
        too_far = False
        j = 0
        while j < len(reverse_primers) and not too_far:
            primer2 = reverse_primers[j]
            if primer2['pos'][0] < primer1['pos'][1]:
                j += 1
                continue
            if primer2['pos'][0] - primer1['pos'][1] > max_overlap:
                too_far = True
                continue

            overlap_start = primer1['pos'][0]
            overlap_end = primer2['pos'][1]
            overlap_length = overlap_end - overlap_start

            if min_overlap <= overlap_length <= max_overlap and primer3.calc_heterodimer_tm(primer1['seq'], primer2['seq']) <= max_heterodimer_tm:
                possible_overlaps.append({
                    'forward': primer1,
                    'reverse': primer2,
                    'overlap_start': overlap_start,
                    'overlap_end': overlap_end,
                    'overlap_length': overlap_length
                })
            j += 1
    
    return possible_overlaps

def get_primer_overlaps_from_seq(seq, allowed_regions, clamp_pattern_forward, clamp_pattern_reverse, clamp_length,
                                 min_length=18, max_length=23, max_poly_x_length=4,
                                 min_gc=0.3, max_gc=0.7, min_tm=57, max_tm=62,
                                 max_hairpin_tm=24, max_homodimer_tm=45, max_3_self_tm=35.0,
                                 min_overlap=60, max_overlap=100, max_heterodimer_tm=45.0,
                                 five_prime_clamp_pattern=None):
    """
    A convenience function that takes a sequence and parameters for primer design and overlap criteria, and returns possible primer overlaps.
    """
    forward_primers = find_forward_primers(seq, allowed_regions, clamp_pattern_forward, clamp_length,
                                           min_length=min_length, max_length=max_length,
                                           max_poly_x_length=max_poly_x_length,
                                           min_gc=min_gc, max_gc=max_gc,
                                           min_tm=min_tm, max_tm=max_tm,
                                           max_hairpin_tm=max_hairpin_tm,
                                           max_homodimer_tm=max_homodimer_tm,
                                           max_3_self_tm=max_3_self_tm,
                                           five_prime_clamp_pattern=five_prime_clamp_pattern)
    
    reverse_primers = find_reverse_primers(seq, allowed_regions, clamp_pattern_reverse, clamp_length,
                                           min_length=min_length, max_length=max_length,
                                           max_poly_x_length=max_poly_x_length,
                                           min_gc=min_gc, max_gc=max_gc,
                                           min_tm=min_tm, max_tm=max_tm,
                                           max_hairpin_tm=max_hairpin_tm,
                                           max_homodimer_tm=max_homodimer_tm,
                                           max_3_self_tm=max_3_self_tm,
                                           five_prime_clamp_pattern=five_prime_clamp_pattern)
    
    possible_overlaps = get_primer_overlaps(forward_primers, reverse_primers,
                                            min_overlap=min_overlap, 
                                            max_overlap=max_overlap,
                                            max_heterodimer_tm=max_heterodimer_tm)
    
    return possible_overlaps
