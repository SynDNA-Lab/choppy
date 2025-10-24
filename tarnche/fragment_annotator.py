from dataclasses import dataclass
from typing import Optional
import re
import networkx as nx
from tqdm import tqdm
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
        overlaps = self.get_possible_overlaps()
        print(f"Constructing overlap graph...")
        self.graph = self.construct_overlap_graph(overlaps)

    def get_possible_overlaps(self) -> list[tuple[int, int]]:
        """
        Generate possible overlap positions based on no-homology regions and boundary motif,
        if provided.

        Returns:
            list[tuple[int, int]]: List of (start, end) tuples for possible overlaps.
        """
        if self.config.motif is None:
            # for free overlaps, we don't care about max_overlap
            # the smaller the overlap, the better
            print("No boundary motif provided, generating free overlaps inside no-homology regions...")
            overlaps = self.get_possible_free_overlaps()
        else:
            print(f"Boundary motif '{self.config.motif}' provided, generating motif-based overlaps...")
            overlaps = self.get_possible_motif_overlaps()

        overlaps.sort(key=lambda x: x[0])
        overlaps = list(dict.fromkeys(overlaps))
        overlaps.insert(0, (0, 0))
        overlaps.append((self.seq_len, self.seq_len))
        return overlaps

    def get_possible_free_overlaps(self):
        overlaps = []
        for region in self.nohom_regions:
            for pos in range(region[0], region[1] - self.config.min_overlap + 1, self.config.min_step):
                overlaps.append((pos, pos + self.config.min_overlap))
        return overlaps

    def get_possible_motif_overlaps(self):
        overlaps = []
        motif_len = len(self.config.motif)
        motif_pattern = re.compile(self.config.motif)
        for region in self.nohom_regions:
            motif_positions = [m.start() for m in motif_pattern.finditer(self.seq, pos=region[0], endpos=region[1])]
            for pos in motif_positions:
                next_pos = next(
                    (
                        p
                        for p in motif_positions
                        if p - pos + motif_len >= self.config.min_overlap
                        and p - pos + motif_len <= self.config.max_overlap
                    ),
                    None,
                )
                if next_pos is not None:
                    overlaps.append((pos, next_pos + motif_len))
        return overlaps

    def construct_overlap_graph(self, overlaps):
        """
        Intial graph construction. Assumes that the provided overlaps are sorted and
        connects them based on min and max fragment size constraints.
        
        Returns:
            nx.DiGraph: Directed graph of overlaps as nodes and valid fragment connections as edges.
        """
        graph = nx.DiGraph()
        graph.add_node((0, 0))
        graph.add_node((self.seq_len, self.seq_len))

        for i in range(len(overlaps)):
            too_far = False
            j = i + 1
            while not too_far and j < len(overlaps):
                if overlaps[j][1] - overlaps[i][0] >= self.config.min_size:
                    if overlaps[j][1] - overlaps[i][0] <= self.config.max_size:
                        graph.add_edge(overlaps[i], overlaps[j], weight=1, constraints="all")
                    else:
                        too_far = True
                j += 1
        
        return graph
    
    def keep_homology_constraint(self, components):
        extra_vertices = []
        for i in range(len(components) - 1):
            gap_start_overlap = max(components[i], key=lambda x: x[1])
            print(gap_start_overlap)
            gap_start = max(gap_start_overlap[0] - self.config.min_size, 0)
            gap_end_overlap = min(components[i+1], key=lambda x: x[0])
            print(gap_end_overlap)
            gap_end = min(gap_end_overlap[1] + self.config.min_size, self.seq_len)

            for region in self.nohom_regions:
                if region[0] < gap_end and region[1] > gap_start:
                    start = max(gap_start, region[0])
                    while start + self.config.min_overlap <= min(region[1], gap_end):
                        extra_vertices.append((start, start + self.config.min_overlap))
                        start += self.config.min_step
        self.add_extra_vertices(extra_vertices, 10, "homology")
    
    def add_extra_vertices(self, extra_vertices, weight, constraints):
        for v in extra_vertices:
            self.graph.add_node(v)
            for u in self.graph.nodes:
                if u[1] < v[0]:
                    if v[1] - u[0] >= self.config.min_size and v[1] - u[0] <= self.config.max_size:
                        self.graph.add_edge(u, v, weight=weight, constraints=constraints)
                elif u[0] > v[1]:
                    if u[1] - v[0] >= self.config.min_size and u[1] - v[0] <= self.config.max_size:
                        self.graph.add_edge(v, u, weight=weight, constraints=constraints)

    def keep_no_constraint(self, components):
        extra_vertices = []
        for i in range(len(components) - 1):
            gap_start_overlap = max(components[i], key=lambda x: x[1])
            gap_start = gap_start_overlap[0]
            gap_end_overlap = min(components[i+1], key=lambda x: x[0])
            gap_end = gap_end_overlap[1]

            # if gap can't fit a fragment, extend it
            if gap_end - gap_start < self.config.min_size:
                gap_start = max(gap_start - self.config.min_size, 0)
                gap_end = min(gap_end + self.config.min_size, self.seq_len)

            step = (self.config.max_size - self.config.min_overlap) // 5
            start = gap_start
            while start + self.config.min_overlap <= gap_end:
                extra_vertices.append((start, start + self.config.min_overlap))
                start += step

        self.add_extra_vertices(extra_vertices, 100, "none")
    
    def fix_graph(self):
        components = list(nx.weakly_connected_components(self.graph))
        components.sort(key=lambda x: min(x, key=lambda y: y[0]))

        print(f"""
            The sequence can't be fully covered with the given constraints.
            Attempting to fix the graph by relaxing constraints...
            """)
        if self.config.motif is not None:
            self.keep_homology_constraint(components)

        fixed = nx.is_weakly_connected(self.graph)
        if not fixed:
            print(""""Failed to fix the graph by keeping homology constraint. 
                    The gap will be closed by unconstrained fragments.
                """)
            self.keep_no_constraint(list(nx.weakly_connected_components(self.graph)))
        else:
            print("Successfully fixed the graph by keeping homology constraint.")

    def get_shortest_path(self):
        if not nx.is_weakly_connected(self.graph):
            self.fix_graph()
        return nx.shortest_path(self.graph, (0, 0), (self.seq_len, self.seq_len), weight="weight")
    
    def get_fragments(self):
        path = self.get_shortest_path()

        fragments = []
        for i in range(len(path) - 1):
            frag_start = path[i][0]
            frag_end = path[i + 1][1]
            overlap_start = path[i + 1][0]
            overlap_end = path[i + 1][1]
            overlap_constraints = self.graph.edges[path[i], path[i + 1]]["constraints"]
            fragments.append((frag_start, frag_end, overlap_start, overlap_end, overlap_constraints))

        return fragments

# end of Fragmentor class

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

def annotate_fragments(record: SeqRecord, config):
    fragmentor = Fragmentor(str(record.seq), extract_no_homology_regions(record), config)
    fragments = fragmentor.get_fragments()

    new_features = list(getattr(record, "features", []))
    frag_records = []
    for idx, (
        frag_start,
        frag_end,
        overlap_start,
        overlap_end,
        overlap_constraints,
    ) in enumerate(fragments, 1):
        frag_feat = SeqFeature(
            FeatureLocation(frag_start, frag_end),
            type="fragment",
            qualifiers={
                "note": f"Fragment_{idx}",
                "overlap_region": f"{overlap_start}-{overlap_end}"
                if overlap_start is not None
                else "None",
                "constraints": overlap_constraints,
            },
        )
        new_features.append(frag_feat)
        frag_seq = record.seq[frag_start:frag_end]
        frag_rec = SeqRecord(
            frag_seq,
            id=f"{record.id}_fragment_{idx}",
            description=f"Fragment {idx}: {frag_start}-{frag_end}, overlap: {overlap_start}-{overlap_end}, constraints: {overlap_constraints}",
        )
        frag_records.append(frag_rec)
    record.features = new_features
    return (record, frag_records)

def fragment_from_file(
    input_file: str,
    output_gb: str,
    output_fasta: str,
    config: FragmentConfig,
):
    input_record = SeqIO.read(input_file, "genbank")
    nohom_regions = extract_no_homology_regions(input_record)
    if not nohom_regions:
        raise RuntimeError("No no-homology regions found in input annotation!")

    (record, frag_records) = annotate_fragments(input_record, config)
    with open(output_gb, "w") as gbo:
        SeqIO.write(record, gbo, "genbank")
    with open(output_fasta, "w") as fao:
        SeqIO.write(frag_records, fao, "fasta")

    print(f"Wrote {len(frag_records)} fragments to {output_gb} and {output_fasta}")
