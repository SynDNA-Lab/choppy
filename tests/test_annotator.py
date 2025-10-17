"""
Tests for the fragment annotator module.
"""
import pytest
import networkx as nx
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from tarnche.fragment_annotator import (
    FragmentConfig,
    extract_no_homology_regions,
    Fragmentor
)


def get_test_data_path(filename):
    """Get the path to a test data file."""
    return Path(__file__).parent / "data" / filename


def load_test_records():
    """Load all test records from the test GenBank file."""
    test_file = get_test_data_path("test_extract_homology_regions.gb")
    records = {}
    for record in SeqIO.parse(test_file, "genbank"):
        records[record.id] = record
    return records


# Load test records once at module level
TEST_RECORDS = load_test_records()

class TestExtractNoHomologyRegions:
    """Tests for extract_no_homology_regions function."""
    
    def test_extract_single_region(self):
        """Test extracting a single non-homology region."""
        record = TEST_RECORDS["single_region"]
        
        regions = extract_no_homology_regions(record)
        
        assert regions == [(5, 10)]
    
    def test_extract_multiple_regions(self):
        """Test extracting multiple non-homology regions."""
        record = TEST_RECORDS["multiple_regions"]
        
        regions = extract_no_homology_regions(record)
        
        assert regions == [(5, 10), (15, 20), (25, 28)]
    
    def test_extract_no_regions(self):
        """Test when there are no homology_free regions."""
        record = TEST_RECORDS["no_regions"]
        
        regions = extract_no_homology_regions(record)
        
        assert regions == []
    
    def test_sort_unordered_regions(self):
        """Test that regions are sorted by start position."""
        record = TEST_RECORDS["unordered_regions"]
        
        regions = extract_no_homology_regions(record)
        
        assert regions == [(5, 10), (12, 18), (20, 25)]
    
    def test_ignore_non_homology_free_features(self):
        """Test that only homology_free features are extracted."""
        record = TEST_RECORDS["mixed_features"]
        
        regions = extract_no_homology_regions(record)
        
        # Only the first misc_feature with homology_free label should be extracted
        assert regions == [(5, 10)]


class TestFragmentorInitialization:
    """Tests for Fragmentor class initialization."""
    
    def test_init_basic(self):
        """Test basic initialization of Fragmentor."""
        seq = "ATGCatgcATGC"
        nohom_regions = [(5, 10)]
        config = FragmentConfig(100, 500, 20, 50)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        assert fragmentor.seq == "ATGCATGCATGC"  # Should be uppercase
        assert fragmentor.seq_len == len(seq)
        assert fragmentor.config == config
        assert fragmentor.nohom_regions == nohom_regions
    
    def test_init_with_multiple_regions(self):
        """Test initialization with multiple no-homology regions."""
        seq = "A" * 1000
        nohom_regions = [(100, 150), (300, 350), (700, 750)]
        config = FragmentConfig(200, 800, 30, 60)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        assert fragmentor.nohom_regions == nohom_regions
        assert fragmentor.seq_len == 1000
    
    def test_init_uppercase_conversion(self):
        """Test that sequence is converted to uppercase."""
        seq = "atgcATGCatgc"
        nohom_regions = []
        config = FragmentConfig(100, 500, 20, 50)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        assert fragmentor.seq == "ATGCATGCATGC"

class TestFragmentorGetPossibleFreeOverlaps:
    """Tests for Fragmentor.get_possible_free_overlaps method."""
    
    def test_free_overlaps_single_region(self):
        """Test generating free overlaps for a single region."""
        seq = "A" * 200
        nohom_regions = [(50, 100)]
        config = FragmentConfig(100, 500, 20, 50, min_step=10)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_free_overlaps()
        
        # Should generate overlaps from position 50 to 100 - 20, step 10
        # Each overlap is (pos, pos + min_overlap)
        expected = [(50, 70), (60, 80), (70, 90), (80, 100)]
        assert overlaps == expected
    
    def test_free_overlaps_multiple_regions(self):
        """Test generating free overlaps for multiple regions."""
        seq = "A" * 300
        nohom_regions = [(50, 80), (150, 180)]
        config = FragmentConfig(100, 500, 21, 50, min_step=10)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_free_overlaps()
        expected = [(50, 71), (150, 171)]
        
        assert overlaps == expected
    
    def test_free_overlaps_small_region(self):
        """Test that small regions don't generate overlaps."""
        seq = "A" * 200
        # Region is only 15bp, smaller than min_overlap (20)
        nohom_regions = [(50, 65)]
        config = FragmentConfig(100, 500, 20, 50, min_step=10)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_free_overlaps()
        
        # Should not generate any overlaps since region is too small
        assert overlaps == []
        
    def test_free_overlaps_exact_min_overlap(self):
        """Test when region size equals min_overlap."""
        seq = "A" * 200
        nohom_regions = [(50, 75)]  # Exactly 20bp
        config = FragmentConfig(100, 500, 20, 50, min_step=10)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_free_overlaps()
        
        # Should generate exactly one overlap at the end
        assert len(overlaps) == 1
        assert overlaps[0][1] - overlaps[0][0] == 20


class TestFragmentorGetPossibleMotifOverlaps:
    """Tests for Fragmentor.get_possible_motif_overlaps method."""
    
    def test_motif_overlaps_found(self):
        """Test generating overlaps based on motif positions."""
        # Create sequence with EcoRI sites (GAATTC)
        seq = "ATGC" + "GAATTC" + "ATGC" * 5 + "GAATTC" + "ATGC" * 5 + "GAATTC" + "ATGC"
        nohom_regions = [(0, len(seq))]
        config = FragmentConfig(100, 500, 20, 40, motif="GAATTC")
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_motif_overlaps()
        
        # Should find pairs of GAATTC positions that meet overlap constraints
        assert len(overlaps) == 2
        # Verify each overlap starts and ends with the motif
        for start, end in overlaps:
            assert seq[start:start+6] == "GAATTC"
            assert seq[end-6:end] == "GAATTC"
    
    def test_motif_overlaps_no_motif_found(self):
        """Test when motif doesn't exist in sequence."""
        seq = "ATGCATGCATGCATGC"
        nohom_regions = [(0, len(seq))]
        config = FragmentConfig(100, 500, 20, 40, motif="GAATTC")
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_motif_overlaps()
        
        assert overlaps == []
    
    def test_motif_overlaps_insufficient_spacing(self):
        """Test when motif occurrences are too close together."""
        # The only possible overlap is 16bp long (including the motifs), but min_overlap is 20
        seq = "ATGC" + "GAATTC" + "ATGC" + "GAATTC" + "ATGC"
        nohom_regions = [(0, len(seq))]
        config = FragmentConfig(100, 500, 20, 40, motif="GAATTC")
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_motif_overlaps()
        
        # Should not find valid overlaps since spacing is insufficient
        assert overlaps == []
    
    def test_motif_overlaps_excessive_spacing(self):
        """Test when motif occurrences are too far apart."""
        # Two GAATTC motifs 60bp apart, but max_overlap is 40 
        seq = "ATGC" + "GAATTC" + "A" * 60 + "GAATTC" + "ATGC"
        nohom_regions = [(0, len(seq))]
        config = FragmentConfig(100, 500, 20, 40, motif="GAATTC")
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_motif_overlaps()
        
        # Should not find valid overlaps since spacing exceeds max_overlap
        assert overlaps == []
    
    def test_motif_overlaps_multiple_regions(self):
        """Test motif overlaps across multiple regions."""
        seq = "GAATTC" + "A" * 20 + "GAATTC" + "T" * 20 + "GAATTC" + "A" * 20 + "GAATTC"
        nohom_regions = [(0, 40), (50, len(seq))]
        config = FragmentConfig(100, 500, 20, 40, motif="GAATTC")
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_motif_overlaps()
        
        # Should find overlaps in both 
        # The overlap that spans across regions is invalid
        assert overlaps == [(0, 32), (52, 84)]
    
    def test_motif_overlaps_exact_boundaries(self):
        """Test motif at exact min and max overlap boundaries."""
        # Create sequence with motifs at exactly 20bp and 40bp apart
        seq = "GAATTC" + "A" * 8 + "GAATTC" + "A" * 7 + "GAATTC" + "A" * 28 + "GAATTC"  + "A" * 29 + "GAATTC"
        nohom_regions = [(0, len(seq))]
        config = FragmentConfig(100, 500, 20, 40, motif="GAATTC")
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_motif_overlaps()
        
        # Should find overlaps at boundary conditions
        for start, end in overlaps:
            assert end - start <= 40
            assert end - start >= 20
        
        assert (0, 20) in overlaps # the smallest overlap that should be found
        assert (27, 67) in overlaps # the largest overlap that should be found


class TestFragmentorGetPossibleOverlaps:
    """Tests for Fragmentor.get_possible_overlaps method."""
    
    def test_get_possible_overlaps_free(self):
        """Test that get_possible_overlaps uses free overlaps when no motif."""
        seq = "A" * 200
        nohom_regions = [(50, 100)]
        config = FragmentConfig(100, 500, 20, 50, min_step=10)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_overlaps()
        
        # Should have start and end sentinels
        assert overlaps[0] == (0, 0)
        assert overlaps[-1] == (len(seq), len(seq))
        # Should have overlaps in between
        assert len(overlaps) == 6  # 4 from region + 2 sentinels
    
    def test_get_possible_overlaps_motif(self):
        """Test that get_possible_overlaps uses motif overlaps when motif provided."""
        seq = "ATGC" + "GAATTC" + "ATGC" * 8 + "GAATTC" + "ATGC"
        nohom_regions = [(0, len(seq))]
        config = FragmentConfig(100, 500, 20, 50, motif="GAATTC")
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_overlaps()
        
        # Should have start and end sentinels
        assert overlaps[0] == (0, 0)
        assert overlaps[-1] == (len(seq), len(seq))
        assert len(overlaps) == 3  # 1 motif overlaps + 2 sentinels
    
    def test_overlaps_sorted(self):
        """Test that overlaps are sorted by start position."""
        seq = "A" * 300
        nohom_regions = [(50, 100), (70, 170), (150, 180)]  # Unordered regions
        config = FragmentConfig(100, 500, 20, 50, min_step=10)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_overlaps()
        
        # Check that overlaps are sorted
        for i in range(len(overlaps) - 1):
            assert overlaps[i][0] <= overlaps[i+1][0]
    
    def test_overlaps_deduplicated(self):
        """Test that duplicate overlaps are removed."""
        seq = "A" * 200
        nohom_regions = [(50, 100), (50, 100)]  # Duplicate regions
        config = FragmentConfig(100, 500, 20, 50, min_step=10)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_overlaps()
        
        # Check no duplicates
        assert len(overlaps) == len(set(overlaps))
class TestFragmentorConstructOverlapGraph:
    """Tests for Fragmentor.construct_overlap_graph method."""
    
    def test_construct_simple_graph(self):
        """Test constructing a simple overlap graph."""
        seq = "A" * 500
        nohom_regions = [(50, 150), (250, 350)]
        config = FragmentConfig(100, 300, 20, 50, min_step=20)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_overlaps()
        graph = fragmentor.construct_overlap_graph(overlaps)
        
        # Graph should be a directed graph
        assert isinstance(graph, nx.DiGraph)
        assert len(graph.nodes) > 0
        assert len(graph.edges) > 0
    
    def test_graph_respects_min_length(self):
        """Test that graph edges constraints."""
        seq = "A" * 500
        nohom_regions = [(50, 110), (99, 130), (159, 250), (200, 300)]
        config = FragmentConfig(100, 300, 20, 50)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        overlaps = fragmentor.get_possible_overlaps()
        graph = fragmentor.construct_overlap_graph(overlaps)
        
        # Check all edges respect min_size and max_size constraints
        for u, v in graph.edges:
            fragment_length = v[1] - u[0]
            assert fragment_length >= 100
            assert fragment_length <= 300
            assert u[0] < v[0]  # Start position must increase
            assert u[1] <= v[1]  # End position should not decrease            

        assert ((0, 0), (80, 100)) in graph.edges # The shortes valid edge
        assert ((200, 220), (500, 500)) in graph.edges # The longest valid edge

class TestFragmentorFixGraph:
    """Tests for Fragmentor.fix_graph and related methods."""
    
    def test_fix_graph_disconnected_graph(self):
        """Test fixing a disconnected graph."""
        # Create a scenario where the graph will be disconnected
        # This happens when there's a gap that can't be covered with the constraints
        seq = "A" * 1000
        # Regions far apart with constraints that can't bridge the gap
        nohom_regions = [(50, 100), (800, 850)]
        config = FragmentConfig(100, 200, 20, 50, min_step=30)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        # Initial graph should be disconnected
        assert not nx.is_weakly_connected(fragmentor.graph)
        
        fragmentor.fix_graph()
        # After fixing, graph should be connected
        assert nx.is_weakly_connected(fragmentor.graph)
        # New edges also respect size constraints
        for u, v in fragmentor.graph.edges:
            fragment_length = v[1] - u[0]
            assert fragment_length >= 100
            assert fragment_length <= 300
            assert u[0] < v[0]  # Start position must increase
            assert u[1] <= v[1]  # End position should not decrease
    
    def test_fix_graph_already_connected(self):
        """Test fix_graph on an already connected graph."""
        seq = "A" * 500
        nohom_regions = [(50, 150), (200, 300)]
        config = FragmentConfig(100, 400, 20, 50, min_step=20)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        # If graph is already connected, fix_graph should not add any edges
        initial_edges = len(fragmentor.graph.edges)
        fragmentor.fix_graph()
        assert nx.is_weakly_connected(fragmentor.graph)
        assert len(fragmentor.graph.edges) == initial_edges
    
    def test_keep_homology_constraint(self):
        """Test that keep_homology_constraint adds vertices in homology regions."""
        seq = "A" * 100 + "GAATTC" + "A" * 20 + "GAATTC" + "A" * 400 + "GAATTC" + "A" * 20 + "GAATTC" + "A" * 100
        # Middle non-homology region has no boundary motifs and can't be bridged
        nohom_regions = [(90, 150), (200, 400), (500, len(seq))]
        config = FragmentConfig(100, 250, 20, 50, min_step=20, motif="GAATTC")
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        initial_nodes = list(fragmentor.graph.nodes)
        
        fragmentor.fix_graph()

        assert nx.is_weakly_connected(fragmentor.graph)
        # Some nodes are added
        new_nodes = set(list(fragmentor.graph.nodes)) - set(initial_nodes)
        assert len(new_nodes) > 0
        # All new nodes are within no-homology regions
        for node in new_nodes:
            print(node)
            in_nohomology = False
            for region in fragmentor.nohom_regions:
                if node[0] >= region[0] and node[1] <= region[1]:
                    in_nohomology = True
                    break
            assert in_nohomology

        # New nodes are at least min_overlap long
        for node in new_nodes:
            assert (node[1] - node[0]) >= fragmentor.config.min_overlap

        # All new edges should have weight=10 and constraints="homology"
        for u, v, data in fragmentor.graph.edges(data=True):
            if u in new_nodes or v in new_nodes:
                assert data.get('constraints') == 'homology'
                assert data['weight'] == 10

        # New edges also respect size constraints
        for u, v in fragmentor.graph.edges:
            fragment_length = v[1] - u[0]
            assert fragment_length >= 100
            assert fragment_length <= 300
            assert u[0] < v[0]  # Start position must increase
            assert u[1] <= v[1]  # End position should not decrease
            
    def test_keep_no_constraint(self):
        """Test that keep_no_constraint adds unconstrained vertices."""
        seq = "A" * 1000
        nohom_regions = [(50, 100), (900, 950)]
        config = FragmentConfig(100, 250, 20, 50, min_step=20)
        
        fragmentor = Fragmentor(seq, nohom_regions, config)
        initial_nodes = list(fragmentor.graph.nodes)
        
        assert not nx.is_weakly_connected(fragmentor.graph)
        fragmentor.fix_graph()

        assert nx.is_weakly_connected(fragmentor.graph)
        # Some nodes are added
        new_nodes = set(list(fragmentor.graph.nodes)) - set(initial_nodes)
        assert len(new_nodes) > 0
        # All new nodes are within no-homology regions

        # New nodes are at least min_overlap long
        for node in new_nodes:
            assert (node[1] - node[0]) >= fragmentor.config.min_overlap

        # All new edges should have weight=100 and constraints="none"
        for u, v, data in fragmentor.graph.edges(data=True):
            if u in new_nodes or v in new_nodes:
                assert data.get('constraints') == 'none'
                assert data['weight'] == 100

    def test_add_extra_vertices(self):
        """Test that add_extra_vertices correctly adds nodes and edges."""
        seq = "A" * 500
        nohom_regions = [(50, 150)]
        config = FragmentConfig(100, 300, 20, 50)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # overlaps = fragmentor.get_possible_overlaps()
        # fragmentor.graph = fragmentor.construct_overlap_graph(overlaps)
        
        # initial_nodes = len(fragmentor.graph.nodes)
        # extra_vertices = [(200, 220), (250, 270)]
        
        # fragmentor.add_extra_vertices(extra_vertices, weight=5, constraints="test")
        
        # Should have added the vertices
        # assert len(fragmentor.graph.nodes) == initial_nodes + 2
        # assert (200, 220) in fragmentor.graph.nodes
        # assert (250, 270) in fragmentor.graph.nodes
        
        # New edges should have the specified weight and constraints
        # for u, v, data in fragmentor.graph.edges(data=True):
        #     if v in extra_vertices or u in extra_vertices:
        #         if data.get('constraints') == 'test':
        #             assert data['weight'] == 5
    
    def test_add_extra_vertices_respects_constraints(self):
        """Test that add_extra_vertices only adds valid edges."""
        seq = "A" * 500
        nohom_regions = [(50, 150)]
        config = FragmentConfig(100, 300, 20, 50)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # overlaps = fragmentor.get_possible_overlaps()
        # fragmentor.graph = fragmentor.construct_overlap_graph(overlaps)
        
        # Add a vertex that would create too-long fragments
        # extra_vertices = [(450, 470)]
        # fragmentor.add_extra_vertices(extra_vertices, weight=5, constraints="test")
        
        # Check that edges respect min_size and max_size
        # for u, v in fragmentor.graph.edges:
        #     fragment_length = v[1] - u[0]
        #     assert fragment_length >= config.min_size
        #     assert fragment_length <= config.max_size
    
    def test_fix_graph_homology_first(self):
        """Test that fix_graph tries homology constraint before no constraint."""
        seq = "A" * 1000
        # Create a gap that CAN be fixed with homology constraint
        nohom_regions = [(100, 400), (600, 900)]
        config = FragmentConfig(150, 600, 30, 60, min_step=30)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # If the graph is disconnected, fix_graph should fix it
        # initial_connected = nx.is_weakly_connected(fragmentor.graph)
        
        # fragmentor.fix_graph()
        
        # Should be connected after fix
        # assert nx.is_weakly_connected(fragmentor.graph)
        
        # Check if homology constraint was used (weight=10 edges exist)
        # has_homology_edges = any(
        #     data.get('constraints') == 'homology' 
        #     for _, _, data in fragmentor.graph.edges(data=True)
        # )
        # assert has_homology_edges
    
    def test_fix_graph_falls_back_to_no_constraint(self):
        """Test that fix_graph falls back to no constraint if homology fails."""
        seq = "A" * 2000
        # Create a gap that CANNOT be fixed with homology constraint alone
        # No homology regions in the gap area
        nohom_regions = [(50, 100), (1900, 1950)]
        config = FragmentConfig(100, 200, 20, 50, min_step=30)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        
        # fragmentor.fix_graph()
        
        # Should be connected after fix
        # assert nx.is_weakly_connected(fragmentor.graph)
        
        # Should have no-constraint edges (weight=100)
        # has_no_constraint_edges = any(
        #     data.get('constraints') == 'none'
        #     for _, _, data in fragmentor.graph.edges(data=True)
        # )
        # assert has_no_constraint_edges


class TestFragmentorShortestPath:
    """Tests for Fragmentor.get_shortest_path method."""
    
    def test_get_shortest_path_simple(self):
        """Test finding shortest path in a simple connected graph."""
        seq = "A" * 500
        nohom_regions = [(50, 150), (250, 350)]
        config = FragmentConfig(100, 300, 20, 50, min_step=30)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # path = fragmentor.get_shortest_path()
        
        # Path should start at (0, 0) and end at (seq_len, seq_len)
        # assert path[0] == (0, 0)
        # assert path[-1] == (500, 500)
        # assert len(path) >= 2
    
    def test_get_shortest_path_disconnected_graph(self):
        """Test shortest path when graph needs fixing."""
        seq = "A" * 1000
        # Create scenario that will need graph fixing
        nohom_regions = [(50, 100), (800, 850)]
        config = FragmentConfig(100, 200, 20, 50, min_step=40)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # Graph might be disconnected initially
        # path = fragmentor.get_shortest_path()
        
        # Should still find a path after fixing
        # assert path[0] == (0, 0)
        # assert path[-1] == (1000, 1000)
        # Path should be continuous (each node connects to next)
        # for i in range(len(path) - 1):
        #     assert (path[i], path[i+1]) in fragmentor.graph.edges
    
    def test_shortest_path_minimizes_weight(self):
        """Test that shortest path minimizes total weight."""
        seq = "A" * 600
        nohom_regions = [(50, 200), (350, 500)]
        config = FragmentConfig(100, 400, 20, 50, min_step=20)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # path = fragmentor.get_shortest_path()
        
        # Calculate path weight
        # path_weight = sum(
        #     fragmentor.graph[path[i]][path[i+1]]['weight']
        #     for i in range(len(path) - 1)
        # )
        
        # Path weight should be reasonable (not excessively high)
        # If we only use "all" constraint edges (weight=1), the weight should be low
        # assert path_weight >= len(path) - 1  # At minimum, 1 per edge
    
    def test_shortest_path_prefers_lower_weight(self):
        """Test that shortest path prefers edges with lower weight."""
        seq = "A" * 500
        nohom_regions = [(50, 150), (250, 350)]
        config = FragmentConfig(100, 300, 20, 50, min_step=20)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # path = fragmentor.get_shortest_path()
        
        # Count edge types in path
        # edge_types = {
        #     'all': 0,
        #     'homology': 0,
        #     'none': 0
        # }
        # for i in range(len(path) - 1):
        #     constraint = fragmentor.graph[path[i]][path[i+1]].get('constraints', 'unknown')
        #     if constraint in edge_types:
        #         edge_types[constraint] += 1
        
        # Should prefer 'all' constraint edges (weight=1) over others
        # if edge_types['all'] > 0:
        #     # If 'all' edges are available, they should dominate
        #     assert edge_types['all'] >= edge_types['homology']
        #     assert edge_types['all'] >= edge_types['none']
    
    def test_shortest_path_covers_entire_sequence(self):
        """Test that the shortest path covers the entire sequence."""
        seq = "A" * 800
        nohom_regions = [(100, 200), (400, 500), (600, 700)]
        config = FragmentConfig(150, 350, 25, 60, min_step=25)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # path = fragmentor.get_shortest_path()
        
        # Verify complete coverage
        # assert path[0] == (0, 0)
        # assert path[-1] == (800, 800)
        
        # Each consecutive pair in the path should overlap or be adjacent
        # for i in range(len(path) - 1):
        #     # Fragment i goes from path[i][0] to path[i+1][1]
        #     # Fragment i+1 starts at path[i+1][0]
        #     # There should be overlap: path[i+1][0] should be between path[i][0] and path[i+1][1]
        #     assert path[i+1][0] <= path[i+1][1]
    
    def test_shortest_path_empty_nohom_regions(self):
        """Test shortest path with no homology regions."""
        seq = "A" * 300
        nohom_regions = []
        config = FragmentConfig(100, 250, 20, 50)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # path = fragmentor.get_shortest_path()
        
        # Should still find a path from start to end
        # assert path[0] == (0, 0)
        # assert path[-1] == (300, 300)
    
    def test_shortest_path_single_fragment(self):
        """Test when entire sequence can be one fragment."""
        seq = "A" * 150
        nohom_regions = [(30, 120)]
        config = FragmentConfig(100, 200, 20, 50)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # path = fragmentor.get_shortest_path()
        
        # Minimum path should have start, one intermediate, and end
        # Or if the entire sequence fits in one fragment: start -> end
        # assert len(path) >= 2
        # assert path[0] == (0, 0)
        # assert path[-1] == (150, 150)
    
    def test_shortest_path_validates_fragment_sizes(self):
        """Test that fragments in the shortest path respect size constraints."""
        seq = "A" * 1000
        nohom_regions = [(100, 300), (500, 700)]
        config = FragmentConfig(150, 400, 30, 70, min_step=30)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # path = fragmentor.get_shortest_path()
        
        # Check each fragment in the path
        # for i in range(len(path) - 1):
        #     fragment_start = path[i][0]
        #     fragment_end = path[i+1][1]
        #     fragment_length = fragment_end - fragment_start
        #     
        #     # Fragment should respect size constraints
        #     assert fragment_length >= config.min_size
        #     assert fragment_length <= config.max_size
    
    def test_shortest_path_with_motif_constraint(self):
        """Test shortest path when using motif-based overlaps."""
        seq = "ATGC" + "GAATTC" + "ATGC" * 30 + "GAATTC" + "ATGC" * 30 + "GAATTC" + "ATGC"
        nohom_regions = [(0, len(seq))]
        config = FragmentConfig(100, 300, 20, 50, motif="GAATTC", min_step=20)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # path = fragmentor.get_shortest_path()
        
        # Should find a path
        # assert path[0] == (0, 0)
        # assert path[-1] == (len(seq), len(seq))
        
        # Verify overlaps align with motif when using intermediate nodes
        # for i in range(1, len(path) - 1):
        #     node = path[i]
        #     # Overlap nodes should be at motif positions
        #     if node[0] > 0 and node[0] < len(seq):
        #         # Check if position aligns with motif
        #         pass  # Implementation depends on exact motif positioning logic


class TestFragmentorEdgeCases:
    """Tests for edge cases in Fragmentor."""
    
    def test_empty_nohom_regions(self):
        """Test with no non-homology regions."""
        seq = "A" * 200
        nohom_regions = []
        config = FragmentConfig(100, 500, 20, 50)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # Should handle empty regions gracefully
    
    def test_very_short_sequence(self):
        """Test with sequence shorter than min fragment size."""
        seq = "ATGC"
        nohom_regions = [(1, 3)]
        config = FragmentConfig(100, 500, 20, 50)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # Should handle short sequences
    
    def test_overlapping_nohom_regions(self):
        """Test with overlapping non-homology regions."""
        seq = "A" * 500
        nohom_regions = [(50, 150), (100, 200)]  # Overlapping
        config = FragmentConfig(100, 300, 20, 50)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # Should handle overlapping regions
    
    def test_adjacent_nohom_regions(self):
        """Test with adjacent non-homology regions."""
        seq = "A" * 500
        nohom_regions = [(50, 150), (150, 250)]  # Adjacent
        config = FragmentConfig(100, 300, 20, 50)
        
        # fragmentor = Fragmentor(seq, nohom_regions, config)
        # Should handle adjacent regions
