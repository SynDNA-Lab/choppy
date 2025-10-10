#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

def extract_no_homology_regions(record):
    """
    Extracts no homology regions (as (start, end) tuples) from the GenBank feature annotations.
    """
    nohom_regions = []
    for feat in getattr(record, 'features', []):
        if feat.type == 'misc_feature':
            if 'label' in feat.qualifiers:
                if feat.qualifiers['label'] == ['homology_free']:
                    start = int(feat.location.start)
                    end = int(feat.location.end)
                    nohom_regions.append((start, end))
    # Sort and merge overlapping/adjacent regions for safety
    nohom_regions.sort()
    merged = []
    for s, e in nohom_regions:
        if not merged or s > merged[-1][1]:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    return [tuple(x) for x in merged]

def region_contains(region, start, end):
    rstart, rend = region
    return start >= rstart and end <= rend

def find_motif_positions(seq, region, motif, is_start=True):
    seq_str = str(seq).upper()
    motif = motif.upper()
    rstart, rend = region
    motif_len = len(motif)
    positions = []
    if is_start:
        for pos in range(rstart, rend - motif_len + 1):
            if seq_str[pos:pos + motif_len] == motif:
                positions.append(pos)
    else:
        for pos in range(rstart + motif_len, rend + 1):
            if seq_str[pos - motif_len:pos] == motif:
                positions.append(pos)
    return positions

def fragment_sequence(seq, nohom_regions, min_size, max_size, min_overlap, max_overlap, motif=None):
    seq_length = len(seq)
    fragments = []
    prev_end = 0
    is_first_fragment = True
    while prev_end < seq_length:
        result = find_next_fragment(
            seq, nohom_regions, prev_end,
            min_size, max_size, min_overlap, max_overlap, motif,
            force_start_at_zero=is_first_fragment
        )
        is_first_fragment = False
        if result is None:
            # Fallback: forcibly advance by min_size (but this will violate constraints)
            frag_start = prev_end
            frag_end = min(prev_end + min_size, seq_length)
            fragments.append((frag_start, frag_end, None, None, None))
            prev_end = frag_end
            continue
        frag_start, frag_end, overlap_start, overlap_end, actual_overlap = result
        fragments.append((frag_start, frag_end, overlap_start, overlap_end, actual_overlap))
        if overlap_start is not None:
            prev_end = overlap_start
        else:
            prev_end = frag_end
    return fragments

def find_next_fragment(seq, nohom_regions, prev_end, min_size, max_size, min_overlap, max_overlap, motif=None, force_start_at_zero=False):
    seq_length = len(seq)
    # If it's the first fragment, always start at 0
    if force_start_at_zero and prev_end == 0:
        start = 0
        # find the largest possible end (max_size) within constraints and motif (if needed)
        for region2 in nohom_regions:
            region2_start = max(region2[0], start + min_size)
            region2_end = min(region2[1], min(start + max_size, seq_length))
            if region2_end <= region2_start:
                continue
            if motif:
                ends = find_motif_positions(seq, (region2_start, region2_end), motif, is_start=False)
            else:
                ends = list(range(region2_start, region2_end + 1))
            for end in sorted(ends, reverse=True):
                frag_len = end - start
                if frag_len < min_size or frag_len > max_size:
                    continue
                # If not last fragment, need overlap in no-homology region and matching motif if set
                if end < seq_length:
                    valid_overlap = False
                    for overlap_region in nohom_regions:
                        ostart = max(overlap_region[0], end - max_overlap)
                        oend = min(overlap_region[1], end)
                        for overlap_len in range(max_overlap, min_overlap - 1, -1):
                            overlap_start = end - overlap_len
                            overlap_end = end
                            if overlap_start < ostart or overlap_end > oend:
                                continue
                            if motif:
                                if str(seq[overlap_start:overlap_start + len(motif)]).upper() != motif.upper():
                                    continue
                            valid_overlap = True
                            return start, end, overlap_start, overlap_end, overlap_len
                    if not valid_overlap:
                        continue
                else:
                    return start, end, None, None, None
        return None  # If nothing is found
    # Normal case: not the very first fragment
    for region in nohom_regions:
        region_start = max(region[0], prev_end)
        region_end = region[1]
        if region_end <= region_start:
            continue
        starts = []
        if motif:
            starts = [p for p in find_motif_positions(seq, (region_start, region_end), motif, is_start=True) if p >= prev_end]
        else:
            starts = list(range(region_start, region_end))
        for start in starts:
            max_frag_end = min(start + max_size, seq_length)
            min_frag_end = start + min_size
            for region2 in nohom_regions:
                if region2[1] <= start + min_size:
                    continue
                region2_start = max(region2[0], min_frag_end)
                region2_end = min(region2[1], max_frag_end)
                if region2_end <= region2_start:
                    continue
                if motif:
                    ends = find_motif_positions(seq, (region2_start, region2_end), motif, is_start=False)
                else:
                    ends = list(range(region2_start, region2_end + 1))
                for end in sorted(ends, reverse=True):
                    frag_len = end - start
                    if frag_len < min_size or frag_len > max_size:
                        continue
                    if end < seq_length:
                        valid_overlap = False
                        for overlap_region in nohom_regions:
                            ostart = max(overlap_region[0], end - max_overlap)
                            oend = min(overlap_region[1], end)
                            for overlap_len in range(max_overlap, min_overlap - 1, -1):
                                overlap_start = end - overlap_len
                                overlap_end = end
                                if overlap_start < ostart or overlap_end > oend:
                                    continue
                                if motif:
                                    if str(seq[overlap_start:overlap_start + len(motif)]).upper() != motif.upper():
                                        continue
                                valid_overlap = True
                                return start, end, overlap_start, overlap_end, overlap_len
                        if not valid_overlap:
                            continue
                    else:
                        return start, end, None, None, None
    return None

def merge_fragments(fragments, max_size):
    """
    Greedily merges consecutive (possibly overlapping) fragments
    as long as the merged fragment does not exceed max_size.
    Preserves the overlap and note fields from the last fragment in the group.
    """
    if not fragments:
        return []
    merged = []
    group_start = fragments[0][0]
    group_end = fragments[0][1]
    group_overlap_start = fragments[0][2]
    group_overlap_end = fragments[0][3]
    group_actual_overlap = fragments[0][4]
    for frag in fragments[1:]:
        next_start, next_end, next_overlap_start, next_overlap_end, next_actual_overlap = frag
        # If merging would exceed max_size, finalize current group
        if next_end - group_start > max_size:
            merged.append((group_start, group_end, group_overlap_start, group_overlap_end, group_actual_overlap))
            group_start = next_start
            group_end = next_end
            group_overlap_start = next_overlap_start
            group_overlap_end = next_overlap_end
            group_actual_overlap = next_actual_overlap
        else:
            # Merge: extend group_end and store latest overlap
            group_end = next_end
            group_overlap_start = next_overlap_start
            group_overlap_end = next_overlap_end
            group_actual_overlap = next_actual_overlap
    # Append last group
    merged.append((group_start, group_end, group_overlap_start, group_overlap_end, group_actual_overlap))
    return merged

def annotate_fragments(record, fragments):
    new_features = list(getattr(record, 'features', []))
    frag_records = []
    for idx, (frag_start, frag_end, overlap_start, overlap_end, actual_overlap) in enumerate(fragments, 1):
        frag_feat = SeqFeature(
            FeatureLocation(frag_start, frag_end),
            type='fragment',
            qualifiers={
                "note": f"Fragment_{idx}",
                "overlap_region": f"{overlap_start}-{overlap_end}" if overlap_start is not None else "None",
                "actual_overlap": str(actual_overlap) if actual_overlap is not None else "None"
            }
        )
        new_features.append(frag_feat)
        frag_seq = record.seq[frag_start:frag_end]
        frag_rec = SeqRecord(
            frag_seq,
            id=f"{record.id}_fragment_{idx}",
            description=f"Fragment {idx}: {frag_start}-{frag_end}, overlap: {overlap_start}-{overlap_end}, actual_overlap: {actual_overlap}"
        )
        frag_records.append(frag_rec)
    record.features = new_features
    return (record, frag_records)

def main():
    parser = argparse.ArgumentParser(description="Fragment an annotated sequence around no homology regions with overlap and boundary constraints.")
    parser.add_argument("input", help="Input GenBank file with no-homology annotated regions (from homology_finder.py)")
    parser.add_argument("--min-size", type=int, required=True, help="Minimum fragment size (bp)")
    parser.add_argument("--max-size", type=int, required=True, help="Maximum fragment size (bp)")
    parser.add_argument("--min-overlap", type=int, required=True, help="Minimum overlap between fragments (bp)")
    parser.add_argument("--max-overlap", type=int, required=True, help="Maximum overlap between fragments (bp)")
    parser.add_argument("--boundary-motif", type=str, default=None, help="Motif that must occur at fragment boundaries (start and end)")
    parser.add_argument("-o", "--output", default="fragmented.gb", help="Output annotated GenBank")
    parser.add_argument("-f", "--fasta", default="fragments.fasta", help="Output multi-fasta file for fragments")
    args = parser.parse_args()

    input_record = SeqIO.read(args.input, "genbank")
    nohom_regions = extract_no_homology_regions(input_record)
    if not nohom_regions:
        raise RuntimeError("No no-homology regions found in input annotation!")

    # Step 1: Fragmentation
    raw_fragments = fragment_sequence(
        input_record.seq, nohom_regions,
        min_size=args.min_size,
        max_size=args.max_size,
        min_overlap=args.min_overlap,
        max_overlap=args.max_overlap,
        motif=args.boundary_motif
    )
    # Step 2: Merge adjacent fragments where possible
    #import pdb; pdb.set_trace()
    merged_fragments = merge_fragments(raw_fragments, args.max_size)

    (record, frag_records) = annotate_fragments(input_record, merged_fragments)
    with open(args.output, "w") as gbo:
        SeqIO.write(record, gbo, "genbank")
    with open(args.fasta, "w") as fao:
        SeqIO.write(frag_records, fao, "fasta")

    print(f"Wrote {len(merged_fragments)} fragments to {args.output} and {args.fasta}")

if __name__ == "__main__":
    main()
