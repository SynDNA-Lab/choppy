#!/usr/bin/env python3

import argparse
import re
import networkx as nx
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord


def extract_no_homology_regions(record: SeqRecord) -> list[tuple[int, int]]:
    """
    Extracts no homology regions (as (start, end) tuples) from the GenBank feature annotations.

    Args:
        record (SeqRecord): A Biopython SeqRecord object with non-homology regions annotated as features.

    Returns:
        list[tuple[int, int]]: A list of (start, end) tuples representing non-homologous regions.
    """
    nohom_regions = []
    for feat in getattr(record, "features", []):
        if feat.type == "misc_feature" and "label" in feat.qualifiers and feat.qualifiers["label"] == ["homology_free"]:
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


def get_possible_overlaps(seq_str, nohom_regions, min_overlap, max_overlap, min_step = 10, motif=None):
    if motif is None:
        # for free overlaps, we don't care about max_overlap
        # the smaller the overlap, the better
        overlaps = get_possible_free_overlaps(nohom_regions, min_overlap, min_step)
    else:
        overlaps = get_possible_motif_overlaps(nohom_regions, min_overlap, max_overlap, motif, seq_str)

    overlaps.sort(key=lambda x: x[0])
    overlaps = list(dict.fromkeys(overlaps))
    overlaps.insert(0, (0, 0))
    overlaps.append((len(seq_str), len(seq_str)))
    return overlaps

def get_possible_free_overlaps(nohom_regions, min_overlap, min_step=10):
    overlaps = []
    for region in nohom_regions:
        for pos in range(region[0] + min_overlap, region[1], min_step):
            overlaps.append((pos - min_overlap, pos))
    return overlaps

def get_possible_motif_overlaps(nohom_regions, min_overlap, max_overlap, motif, seq_str):
    overlaps = []
    motif_len = len(motif)
    for region in nohom_regions:
        motif_positions = [m.start() for m in re.finditer(motif, seq_str, pos=region[0], endpos=region[1])]
        for pos in motif_positions:
            next_pos = next(
                (
                    p
                    for p in motif_positions
                    if p - pos + motif_len - 1 >= min_overlap
                    and p - pos + motif_len - 1 <= max_overlap
                ),
                None,
            )
            if next_pos is not None:
                overlaps.append((pos, next_pos))
    return overlaps

def construct_overlap_graph(overlaps, min_length, max_length):
    graph = nx.DiGraph()
    for i in range(len(overlaps)):
        too_far = False
        j = i + 1
        while not too_far and j < len(overlaps):
            if overlaps[j][1] - overlaps[i][0] >= min_length:
                if overlaps[j][1] - overlaps[i][0] <= max_length:
                    graph.add_edge(overlaps[i], overlaps[j])
                else:
                    too_far = True
            j += 1
    return graph

def fix_graph(graph, nohom_regions, min_length, max_length, min_overlap):
    components = list(nx.weaksy_connected_components(graph))
    comp_span = []
    for comp in components:
        min_node = min(comp, key=lambda x: x[0])
        max_node = max(comp, key=lambda x: x[1])
        comp_span.append((min_node, max_node))
    comp_span.sort(key=lambda x: x[0])
    for i in range(len(comp_span) - 1):
        print(f"""Warning: Gap from position {comp_span[i][1][1]} to {comp_span[i+1][0][0]} 
                cannot be covered by a valid fragment. Gap size: {comp_span[i+1][0][0] - comp_span[i][1][1]}.
                Attempting to fix it by relaxing constraints...
              """)
        success = keep_homology_constraint(graph, (comp_span[i][1][1], comp_span[i+1][0][0]), nohom_regions)
        if not success:
            print(""""Failed to fix the graph by keeping homology constraint. 
                    The gap will be closed by unconstrained fragments.
                """)
            keep_no_constraint(graph, (comp_span[i][1], comp_span[i+1][0]), min_length, max_length, min_overlap)
        else:
            print("Successfully fixed the graph by keeping homology constraint.")
    return graph

def keep_homology_constraint(graph, gap, nohom_regions):
    return True

def keep_no_constraint(graph, gap, min_length, max_length, min_overlap):
    gap_start, gap_end = gap
    graph.add_edge((gap_start, gap_start), (gap_end, gap_end))
    return True

def get_shortest_path(graph):
    if not nx.weakly_connected(graph):
        graph = fix_graph(graph)
    return nx.shortest_path(graph, (0, 0), (len(seq), len(seq)))


def region_contains(region, start, end):
    rstart, rend = region
    return start >= rstart and end <= rend


def annotate_fragments(record, fragments):
    new_features = list(getattr(record, "features", []))
    frag_records = []
    for idx, (
        frag_start,
        frag_end,
        overlap_start,
        overlap_end,
        actual_overlap,
    ) in enumerate(fragments, 1):
        frag_feat = SeqFeature(
            FeatureLocation(frag_start, frag_end),
            type="fragment",
            qualifiers={
                "note": f"Fragment_{idx}",
                "overlap_region": f"{overlap_start}-{overlap_end}"
                if overlap_start is not None
                else "None",
                "actual_overlap": str(actual_overlap)
                if actual_overlap is not None
                else "None",
            },
        )
        new_features.append(frag_feat)
        frag_seq = record.seq[frag_start:frag_end]
        frag_rec = SeqRecord(
            frag_seq,
            id=f"{record.id}_fragment_{idx}",
            description=f"Fragment {idx}: {frag_start}-{frag_end}, overlap: {overlap_start}-{overlap_end}, actual_overlap: {actual_overlap}",
        )
        frag_records.append(frag_rec)
    record.features = new_features
    return (record, frag_records)


def main():
    parser = argparse.ArgumentParser(
        description="Fragment an annotated sequence around no homology regions with overlap and boundary constraints."
    )
    parser.add_argument(
        "input",
        help="Input GenBank file with no-homology annotated regions (from homology_finder.py)",
    )
    parser.add_argument(
        "--min-size", type=int, required=True, help="Minimum fragment size (bp)"
    )
    parser.add_argument(
        "--max-size", type=int, required=True, help="Maximum fragment size (bp)"
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        required=True,
        help="Minimum overlap between fragments (bp)",
    )
    parser.add_argument(
        "--max-overlap",
        type=int,
        required=True,
        help="Maximum overlap between fragments (bp)",
    )
    parser.add_argument(
        "--boundary-motif",
        type=str,
        default=None,
        help="Motif that must occur at fragment boundaries (start and end)",
    )
    parser.add_argument(
        "-o", "--output", default="fragmented.gb", help="Output annotated GenBank"
    )
    parser.add_argument(
        "-f",
        "--fasta",
        default="fragments.fasta",
        help="Output multi-fasta file for fragments",
    )
    args = parser.parse_args()

    input_record = SeqIO.read(args.input, "genbank")
    nohom_regions = extract_no_homology_regions(input_record)
    if not nohom_regions:
        raise RuntimeError("No no-homology regions found in input annotation!")

    # Step 1: Fragmentation
    raw_fragments = fragment_sequence(
        input_record.seq,
        nohom_regions,
        min_size=args.min_size,
        max_size=args.max_size,
        min_overlap=args.min_overlap,
        max_overlap=args.max_overlap,
        motif=args.boundary_motif,
    )
    # Step 2: Merge adjacent fragments where possible
    merged_fragments = merge_fragments(raw_fragments, args.max_size)

    (record, frag_records) = annotate_fragments(input_record, merged_fragments)
    with open(args.output, "w") as gbo:
        SeqIO.write(record, gbo, "genbank")
    with open(args.fasta, "w") as fao:
        SeqIO.write(frag_records, fao, "fasta")

    print(f"Wrote {len(merged_fragments)} fragments to {args.output} and {args.fasta}")


if __name__ == "__main__":
    main()
