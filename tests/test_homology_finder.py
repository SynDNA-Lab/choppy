"""
Tests for the homology finder module.
"""
import pytest
from pathlib import Path
from tarnche.homology_finder import (_parse_sequence_files, _create_kmer_trie)
import marisa_trie as mt

def get_test_data_path(filename):
    """Get the path to a test data file."""
    return Path(__file__).parent / "data" / filename


def test_read_fasta_file():
    """Test reading a FASTA file."""
    fasta_path = get_test_data_path("test_sequence.fa")
    sequences = _parse_sequence_files(str(fasta_path))
    
    assert len(sequences) == 1, "Should have exactly one sequence"
    assert len(sequences[0].seq) == 10, "Sequence should be 10 bases long"
    assert sequences[0].id == "test_seq_1", "Sequence ID should match"


def test_read_genbank_file():
    """Test reading a GenBank file."""
    gb_path = get_test_data_path("test_sequences.gb")
    sequences = _parse_sequence_files(str(gb_path))
    
    assert len(sequences) == 2, "Should have exactly two sequences"
    assert len(sequences[0].seq) == 5, "First sequence should be 5 bases long"
    assert len(sequences[1].seq) == 16, "Second sequence should be 9 bases long"
    assert sequences[0].id == "test_seq_1", "First sequence ID should match"
    assert sequences[1].id == "test_seq_2", "Second sequence ID should match"


def test_read_multiple_files():
    """Test reading both FASTA and GenBank files together."""
    fasta_path = get_test_data_path("test_sequence.fa")
    gb_path = get_test_data_path("test_sequences.gb")
    
    sequences = _parse_sequence_files([str(fasta_path), str(gb_path)])
    
    assert len(sequences) == 3, "Should have three sequences total"
    
    # Check lengths
    lengths = [len(seq.seq) for seq in sequences]
    expected_lengths = [10, 5, 16]  # FASTA (10), GB seq1 (5), GB seq2 (16)
    assert lengths == expected_lengths, f"Expected lengths {expected_lengths}, got {lengths}"


def test_read_csv_file_error():
    """Test that reading a CSV file raises an appropriate error."""
    csv_path = get_test_data_path("test_file.csv")
    
    with pytest.raises(Exception):
        _parse_sequence_files(str(csv_path))


def test_read_nonexistent_file():
    """Test that reading a non-existent file raises an error."""
    with pytest.raises(FileNotFoundError):
        _parse_sequence_files("nonexistent_file.fa")


def test_kmer_trie():
    """Test k-mer trie creation."""
    gb_path = get_test_data_path("test_sequences.gb")
    sequences = _parse_sequence_files(str(gb_path))
    
    # Create k-mer trie for the first sequence with k=3, bg=True
    trie_bg = _create_kmer_trie(sequences[0], kmer_size=4, bg=True)
    assert isinstance(trie_bg, mt.Trie), "Should return a marisa_trie.Trie object"
    assert len(trie_bg) == 9, "Trie should contain 9 unique k-mers (circular sequence)"
    
    # Create k-mer trie for the second sequence with k=4, bg=False
    trie_query = _create_kmer_trie(sequences[1], kmer_size=4, bg=False)
    assert len(trie_query) == 9, "Trie should contain 9 k-mers (query sequence)"

