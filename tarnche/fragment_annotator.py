from dataclasses import dataclass
from typing import Optional
import re
import networkx as nx
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

@dataclass
class FragmentConfig:
    min_size: int
    max_size: int
    min_overlap: int
    max_overlap: int
    motif: Optional[str] = None
    min_step: int = 10

    def __post_init__(self):
        """Validate parameters."""
        if self.min_size > self.max_size:
            raise ValueError("min_size must be <= max_size")
        if self.min_overlap > self.max_overlap:
            raise ValueError("min_overlap must be <= max_overlap")
        # TO DO: Add more checks
class Fragmentor:
    def __init__(self, seq: str, nohom_regions: tuple[int, int], config: FragmentConfig):
        self.seq = seq.upper()
        self.seq_len = len(seq)
        self.config = config
        self.nohom_regions = nohom_regions

        print(f"Calculating possible overlaps...")
        overlaps = self.get_possible_overlaps(self)
        print(f"Constructing overlap graph...")
        graph = self.construct_overlap_graph(self, overlaps)

    def get_possible_overlaps(self):
        if self.config.motif is None:
            # for free overlaps, we don't care about max_overlap
            # the smaller the overlap, the better
            overlaps = self.get_possible_free_overlaps(self)
        else:
            overlaps = self.get_possible_motif_overlaps(self)

        overlaps.sort(key=lambda x: x[0])
        overlaps = list(dict.fromkeys(overlaps))
        overlaps.insert(0, (0, 0))
        overlaps.append((self.seq_len, self.seq_len))
        return overlaps

    def get_possible_free_overlaps(self):
        overlaps = []
        for region in self.nohom_regions:
            for pos in range(region[0] + self.min_overlap, region[1], self.min_step):
                overlaps.append((pos - self.min_overlap, pos))
        return overlaps

    def get_possible_motif_overlaps(self):
        overlaps = []
        motif_len = len(self.motif)
        for region in self.nohom_regions:
            motif_positions = [m.start() for m in re.finditer(self.motif, self.seq, pos=region[0], endpos=region[1])]
            for pos in motif_positions:
                next_pos = next(
                    (
                        p
                        for p in motif_positions
                        if p - pos + motif_len - 1 >= self.min_overlap
                        and p - pos + motif_len - 1 <= self.max_overlap
                    ),
                    None,
                )
                if next_pos is not None:
                    overlaps.append((pos, next_pos))
        return overlaps

    def construct_overlap_graph(self, overlaps):
        graph = nx.DiGraph()
        for i in range(len(overlaps)):
            too_far = False
            j = i + 1
            while not too_far and j < len(overlaps):
                if overlaps[j][1] - overlaps[i][0] >= self.min_length:
                    if overlaps[j][1] - overlaps[i][0] <= self.max_length:
                        graph.add_edge(overlaps[i], overlaps[j], weight=1, constraints="all")
                    else:
                        too_far = True
                j += 1
        return graph
    
    def keep_homology_constraint(self, components):
        extra_vertices = []
        for i in range(len(components) - 1):
            gap_start = max(components[i], key=lambda x: x[1])
            gap_start[0] -= self.min_length
            gap_start[0] = max(gap_start[0], 0)
            gap_end = min(components[i+1], key=lambda x: x[0])
            gap_end[1] += self.min_length
            gap_end[1] = min(gap_end[1], self.seq_len)

            for region in self.nohom_regions:
                if region[0] < gap_end[1] and region[1] > gap_start[0]:
                    start = max(gap_start[1], region[0])
                    while start + self.min_overlap <= min(region[1], gap_end[0]):
                        extra_vertices.append((start - self.min_overlap, start))
                        start += self.min_step
        self.add_extra_vertices(self, extra_vertices, 10, "homology")
    
    def add_extra_vertices(self, extra_vertices, weight, constraints):
        for v in extra_vertices:
            self.graph.add_node(v)
            for u in self.graph.nodes:
                if u[1] < v[0]:
                    if v[1] - u[0] >= self.min_length and v[1] - u[0] <= self.max_length:
                        self.graph.add_edge(u, v, weight=weight, constraints=constraints)
                elif u[0] > v[1]:
                    if u[1] - v[0] >= self.min_length and u[1] - v[0] <= self.max_length:
                        self.graph.add_edge(v, u, weight=weight, constraints=constraints)

    def keep_no_constraint(self, components):
        extra_vertices = []
        for i in range(len(components) - 1):
            gap_start = max(components[i], key=lambda x: x[1])
            gap_start[0] -= self.min_length
            gap_start[0] = max(gap_start[0], 0)
            gap_end = min(components[i+1], key=lambda x: x[0])
            gap_end[1] += self.min_length
            gap_end[1] = min(gap_end[1], self.seq_len)

            start = gap_start[0]
            while start + self.min_overlap <= gap_end[1]:
                extra_vertices.append((start - self.min_overlap, start))
                start += self.max_length

        self.add_extra_vertices(self, extra_vertices, 100, "none")
    
    def fix_graph(self):
        components = list(nx.weakly_connected_components(self.graph))
        components.sort(key=lambda x: min(x, key=lambda y: y[0]))

        print(f"""
            The sequence can't be fully covered with the given constraints.
            Attempting to fix the graph by relaxing constraints...
            """)
        
        self.keep_homology_constraint(self, components)

        fixed = self.graph.is_weakly_connected()
        if not fixed:
            print(""""Failed to fix the graph by keeping homology constraint. 
                    The gap will be closed by unconstrained fragments.
                """)
            self.keep_no_constraint(self, self.graph.weakly_connected_components())
        else:
            print("Successfully fixed the graph by keeping homology constraint.")

    def get_shortest_path(self):
        if not nx.is_weakly_connected(self.graph):
            graph = self.fix_graph()
        return nx.shortest_path(graph, (0, 0), (self.seq_len, self.seq_len), weight="weight")

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
    nohom_regions.sort()
    return nohom_regions

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


# def main():
#     parser = argparse.ArgumentParser(
#         description="Fragment an annotated sequence around no homology regions with overlap and boundary constraints."
#     )
#     parser.add_argument(
#         "input",
#         help="Input GenBank file with no-homology annotated regions (from homology_finder.py)",
#     )
#     parser.add_argument(
#         "--min-size", type=int, required=True, help="Minimum fragment size (bp)"
#     )
#     parser.add_argument(
#         "--max-size", type=int, required=True, help="Maximum fragment size (bp)"
#     )
#     parser.add_argument(
#         "--min-overlap",
#         type=int,
#         required=True,
#         help="Minimum overlap between fragments (bp)",
#     )
#     parser.add_argument(
#         "--max-overlap",
#         type=int,
#         required=True,
#         help="Maximum overlap between fragments (bp)",
#     )
#     parser.add_argument(
#         "--boundary-motif",
#         type=str,
#         default=None,
#         help="Motif that must occur at fragment boundaries (start and end)",
#     )
#     parser.add_argument(
#         "-o", "--output", default="fragmented.gb", help="Output annotated GenBank"
#     )
#     parser.add_argument(
#         "-f",
#         "--fasta",
#         default="fragments.fasta",
#         help="Output multi-fasta file for fragments",
#     )
#     args = parser.parse_args()

#     input_record = SeqIO.read(args.input, "genbank")
#     nohom_regions = extract_no_homology_regions(input_record)
#     if not nohom_regions:
#         raise RuntimeError("No no-homology regions found in input annotation!")

#     # Step 1: Fragmentation
#     raw_fragments = fragment_sequence(
#         input_record.seq,
#         nohom_regions,
#         min_size=args.min_size,
#         max_size=args.max_size,
#         min_overlap=args.min_overlap,
#         max_overlap=args.max_overlap,
#         motif=args.boundary_motif,
#     )
#     # Step 2: Merge adjacent fragments where possible
#     merged_fragments = merge_fragments(raw_fragments, args.max_size)

#     (record, frag_records) = annotate_fragments(input_record, merged_fragments)
#     with open(args.output, "w") as gbo:
#         SeqIO.write(record, gbo, "genbank")
#     with open(args.fasta, "w") as fao:
#         SeqIO.write(frag_records, fao, "fasta")

#     print(f"Wrote {len(merged_fragments)} fragments to {args.output} and {args.fasta}")