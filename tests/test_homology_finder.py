"""
Tests for the homology finder module.
"""
import pytest
from pathlib import Path
from tarnche.homology_finder import (
    parse_sequence_files, 
    create_kmer_trie,
    merge_tries,
    process_background_sequences,
    process_query_sequences,
    find_non_homologous_regions
)
import marisa_trie as mt
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def get_test_data_path(filename):
    """Get the path to a test data file."""
    return Path(__file__).parent / "data" / filename


def test_read_fasta_file():
    """Test reading a FASTA file."""
    fasta_path = get_test_data_path("test_sequence.fa")
    sequences = parse_sequence_files(str(fasta_path))
    
    assert len(sequences) == 1, "Should have exactly one sequence"
    assert len(sequences[0].seq) == 10, "Sequence should be 10 bases long"
    assert sequences[0].id == "test_seq_1", "Sequence ID should match"


def test_read_genbank_file():
    """Test reading a GenBank file."""
    gb_path = get_test_data_path("test_sequences.gb")
    sequences = parse_sequence_files(str(gb_path))
    
    assert len(sequences) == 2, "Should have exactly two sequences"
    assert len(sequences[0].seq) == 5, "First sequence should be 5 bases long"
    assert len(sequences[1].seq) == 16, "Second sequence should be 9 bases long"
    assert sequences[0].id == "test_seq_1", "First sequence ID should match"
    assert sequences[1].id == "test_seq_2", "Second sequence ID should match"


def test_read_multiple_files():
    """Test reading both FASTA and GenBank files together."""
    fasta_path = get_test_data_path("test_sequence.fa")
    gb_path = get_test_data_path("test_sequences.gb")
    
    sequences = parse_sequence_files([str(fasta_path), str(gb_path)])
    
    assert len(sequences) == 3, "Should have three sequences total"
    
    # Check lengths
    lengths = [len(seq.seq) for seq in sequences]
    expected_lengths = [10, 5, 16]  # FASTA (10), GB seq1 (5), GB seq2 (16)
    assert lengths == expected_lengths, f"Expected lengths {expected_lengths}, got {lengths}"


def test_read_csv_file_error():
    """Test that reading a CSV file raises an appropriate error."""
    csv_path = get_test_data_path("test_file.csv")
    
    with pytest.raises(Exception):
        parse_sequence_files(str(csv_path))


def test_read_nonexistent_file():
    """Test that reading a non-existent file raises an error."""
    with pytest.raises(FileNotFoundError):
        parse_sequence_files("nonexistent_file.fa")


def test_kmer_trie():
    """Test k-mer trie creation."""
    gb_path = get_test_data_path("test_sequences.gb")
    sequences = parse_sequence_files(str(gb_path))
    
    # Create k-mer trie for the first sequence with k=3, bg=True
    trie_bg = create_kmer_trie(sequences[0], kmer_size=4, bg=True)
    assert isinstance(trie_bg, mt.Trie), "Should return a marisa_trie.Trie object"
    assert len(trie_bg) == 9, "Trie should contain 9 unique k-mers (circular sequence)"
    
    # Create k-mer trie for the second sequence with k=4, bg=False
    trie_query = create_kmer_trie(sequences[1], kmer_size=4, bg=False)
    assert len(trie_query) == 9, "Trie should contain 9 k-mers (query sequence)"

def test_merge_tries():
    """Test merging of two k-mer tries."""
    gb_path = get_test_data_path("test_sequences.gb")
    sequences = parse_sequence_files(str(gb_path))
    
    trie1 = create_kmer_trie(sequences[0], kmer_size=4, bg=False)
    trie2 = create_kmer_trie(sequences[1], kmer_size=4, bg=False)

    merged_trie = merge_tries([trie1, trie2])
    assert isinstance(merged_trie, mt.Trie), "Merged result should be a marisa_trie.Trie object"
    assert len(merged_trie) == 10, "Merged trie should have 10 unique k-mers"

def test_process_background_sequences():
    """Test processing of background sequences."""
    gb_path = get_test_data_path("test_sequences.gb")
    sequences = parse_sequence_files(str(gb_path))
    
    bg_trie = process_background_sequences(sequences, kmer_size=4)
    assert isinstance(bg_trie, mt.Trie), "Should return a marisa_trie.Trie object"
    assert len(bg_trie) == 22, "Background trie should contain 22 unique k-mers"

    empty_trie = process_background_sequences([], kmer_size=4)
    assert len(empty_trie) == 0, "Empty input should return an empty trie"

def test_process_query_sequences():
    """Test processing of query sequences."""
    gb_path = get_test_data_path("test_sequences.gb")
    sequences = parse_sequence_files(str(gb_path))
    
    query_tries = process_query_sequences(sequences, kmer_size=4)
    assert isinstance(query_tries, dict), "Should return a dict of tries"
    assert len(query_tries) == 2, "dict contains 2 tries"
    assert all(isinstance(trie, mt.Trie) for trie in query_tries.values()), "All values should be marisa_trie.Trie objects"
    assert query_tries.keys() == {sequences[0].id, sequences[1].id}, "Keys should match sequence IDs"
    assert len(query_tries[sequences[0].id]) == 1, "First trie should have 1 k-mers"
    assert len(query_tries[sequences[1].id]) == 9, "Second trie should have 9 k-mers"


def test_find_non_homologous_regions():
    """Test finding non-homologous regions."""
    seq = SeqRecord(Seq("ATGCATGCAAGGCCTT"), id="test", description="test sequence")
    query_trie = create_kmer_trie(seq, kmer_size=4, bg=False)
    bg_trie = mt.Trie()
    regions = find_non_homologous_regions(seq, query_trie, bg_trie, kmer_size=4, threshold=5)
    assert len(regions) == 1, "Should find one non-homologous region"
    assert regions[0] == (6, 11), "Non-homologous region should be from index 6 to 11"

    query_seq = SeqRecord(Seq("ATATAACCCCCTTTTTGGGGCCCGCG"), id="test2", description="test sequence 2")
    bg_seq = SeqRecord(Seq("TATAACCCCCTTTTTGGGG"), id="bg", description="background sequence")
    query_trie = create_kmer_trie(query_seq, kmer_size=5, bg=False)
    bg_trie = create_kmer_trie(bg_seq, kmer_size=5, bg=True)
    regions = find_non_homologous_regions(query_seq, query_trie, bg_trie, kmer_size=5, threshold=5)
    assert len(regions) == 3, "Should find two non-homologous regions"
    assert regions[0] == (0, 5), "First non-homologous region should be from index 0 to 5"
    assert regions[1] == (16, 21), "Second non-homologous region should be from index 16 to 21"
    assert regions[2] == (19, 26), "Third non-homologous region should be from index 19 to 26"
    regions = find_non_homologous_regions(query_seq, query_trie, bg_trie, kmer_size=5, threshold=6)
    assert len(regions) == 1, "Only one region meets the higher threshold"
    assert regions[0] == (19, 26), "Only the last region meets the higher threshold"
