from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from tqdm import tqdm
from typing import List, Dict, Union
import marisa_trie as mt


def parse_sequence_files(file_paths: Union[str, List[str]]) -> List[SeqRecord]:
    """
    Parse one or multiple sequence files (GenBank or FASTA) and return list of SeqRecord objects.

    Args:
        file_paths (Union[str, List[str]]): Path or list of paths to sequence files

    Returns:
        List[SeqRecord]: List of parsed sequence records
    """
    if isinstance(file_paths, str):
        file_paths = [file_paths]

    records = []
    for file_path in file_paths:
        file_ext = Path(file_path).suffix.lower()
        if file_ext not in [".gb", ".gbk", ".genbank", ".fa"]:
            raise Exception(f"Filetype of file '{file_path}' is not supported.")

        format_type = "genbank" if file_ext in [".gb", ".gbk", ".genbank"] else "fasta"
        # Parse all sequences in the file
        records.extend(list(SeqIO.parse(file_path, format_type)))

    return records

# TO DO: add option to subtract background
def create_kmer_trie(sequence: SeqRecord, kmer_size: int, bg: bool) -> mt.Trie:
    """
    Create a k-mer trie from a sequence, considering both strands.

    Args:
        sequence (SeqRecord): Input sequence record
        kmer_size (int): Size of k-mers
        bg (bool): If True, create a set of all unique k-mers (for background sequences); 
            if False, track k-mers appearing more than once (for query sequences)

    Returns:
        mt.Trie: Trie containing the k-mers (marisa_trie)
    """

    sequence_str = str(sequence.seq)
    sequence_rc_str = str(sequence.seq.reverse_complement())

    if sequence.annotations.get("topology", "").lower() == "circular":
        sequence_str = sequence_str + sequence_str[0:kmer_size - 1]
        sequence_rc_str = sequence_rc_str + sequence_rc_str[0:kmer_size - 1]

    if bg:
        kmers = set()
        for i in tqdm(range(len(sequence_str) - kmer_size + 1), desc="Processing sequence"):
            kmers.add(sequence_str[i : i + kmer_size])
            kmers.add(sequence_rc_str[i : i + kmer_size])
    else:
        kmers_dict = {}
        for i in tqdm(range(len(sequence_str) - kmer_size + 1), desc="Processing sequence"):
            for kmer in (sequence_str[i : i + kmer_size], sequence_rc_str[i : i + kmer_size]):
                if kmers_dict.get(kmer) is None:
                    kmers_dict[kmer] = False
                elif kmers_dict[kmer] is False:
                    kmers_dict[kmer] = True
        
        kmers = set(kmer for kmer, val in kmers_dict.items() if val)

    print(f"Found {len(kmers)} unique k-mers of size {kmer_size} in sequence {sequence.id}")
    print("Constructing a trie...")
    return mt.Trie(kmers)

def merge_tries(tries: List[mt.Trie]) -> mt.Trie:
    """
    Merge multiple marisa_tries into a single trie.

    Args:
        tries (List[mt.Trie]): List of marisa_trie.Trie objects to merge

    Returns:
        mt.Trie: Merged marisa_trie.Trie object
    """
    if len(tries) == 0:
        return mt.Trie()
    if len(tries) == 1:
        return tries[0] 
    all_kmers = set()
    for trie in tries:
        all_kmers.update(trie.keys())
    
    return mt.Trie(all_kmers)

def process_background_sequences(
    background_sequences: List[SeqRecord], kmer_size: int, merge: bool = True
) -> mt.Trie:
    """
    Process background sequences to create a combined k-mer trie.

    Args:
        background_sequences (List[SeqRecord]): List of background sequences
        kmer_size (int): Size of k-mers

    Returns:
        mt.Trie: Combined k-mer trie from all background sequences
    """
    tries = []
    for bg_seq in background_sequences:
        trie = create_kmer_trie(bg_seq, kmer_size, bg=True)
        tries.append(trie)

    if merge:
        if tries:
            combined_trie = merge_tries(tries)
        else:
            combined_trie = mt.Trie()
        return combined_trie
    else:
        if len(tries) == 0:
            tries.append(mt.Trie())
        return tries
  

def process_query_sequences(
    query_sequences: List[SeqRecord], kmer_size: int
) -> Dict[str, mt.Trie]:
    """
    Process query sequences to create k-mer tries for each sequence.

    Args:
        query_sequences (List[SeqRecord]): List of query sequences
        kmer_size (int): Size of k-mers 
    
    Returns:
        Dict[str, mt.Trie]: Dictionary of k-mer tries for each query sequence
    """
    results = {}
    for query_seq in query_sequences:
        trie = create_kmer_trie(query_seq, kmer_size, bg=False)
        results[query_seq.id] = trie

    return results

def check_kmer_in_tries(tries: List[mt.Trie], kmer: str) -> bool:
    """
    Check if a k-mer exists in any of the provided tries.

    Args:
        tries (List[mt.Trie]): List of k-mer tries  
        kmer (str): K-mer to check
    Returns:
        bool: True if k-mer is found in any trie, False otherwise
    """
    for trie in tries:
        if kmer in trie:
            return True
    return False

def find_non_homologous_regions(
    sequence: SeqRecord, query_trie: mt.Trie, 
    bg_tries: Union[mt.Trie, List[mt.Trie]], 
    kmer_size: int, threshold: int
) -> list:
    """
    Find non-homologous regions in a sequence based on k-mer tries.

    Args:
        sequence (SeqRecord): Input query sequence
        query_trie (mt.Trie): K-mer trie for the query sequence
        bg_tries (mt.Trie or List[mt.Treis]): K-mer trie or tries for the background sequences
        kmer_size (int): Size of k-mers
        threshold (int): Minimum size of non-homologous regions to report

    Returns:
        List[(int, int)]: List of tuples representing start and end positions of non-homologous regions
    """
    sequence_str = str(sequence.seq)

    if isinstance(bg_tries, mt.Trie):
        bg_tries = [bg_tries]

    if sequence.annotations.get("topology", "").lower() == "circular":
        sequence_str = sequence_str + sequence_str[0:kmer_size - 1]

    regions = []
    start = 0
    opened = False
    for i in tqdm(range(len(sequence_str) - kmer_size + 1), desc="Finding non-homologous regions"):
        kmer = sequence_str[i : i + kmer_size]
        if not check_kmer_in_tries(bg_tries, kmer) and kmer not in query_trie:
            if not opened:
                start = i
                opened = True
        else:
            if opened and i - start + kmer_size - 1 >= threshold:
                regions.append((start, i + kmer_size - 1)) 
            opened = False
    if opened and len(sequence_str) - start >= threshold:
        regions.append((start, len(sequence_str)))

    return regions    

def create_annotated_record(
    sequence: SeqRecord, regions: list
) -> SeqRecord:
    """
    Create an annotated sequence record with non-homologous regions.

    Args:
        sequence (SeqRecord): Input sequence record
        regions (dict): Dictionary of non-homologous regions

    Returns:
        SeqRecord: Annotated sequence record
    """
    features = getattr(sequence, "features", [])

    # Add features for non-homologous regions
    for i, (start, end) in enumerate(regions):
        feature = SeqFeature(
            FeatureLocation(start, end),
            type="misc_feature",
            qualifiers={
                "note": f"Non-homologous region {i + 1}",
                "label": f"homology_free",
            },
        )
        features.append(feature)

    # Create new record
    record = SeqRecord(
        sequence.seq,
        id=f"{sequence.id}_annotated",
        name=sequence.name,
        description="Sequence annotated with non-homologous regions",
        features=features,
        annotations={"molecule_type": "DNA"},
    )

    return record

def file_process_homology(query_path, background_path, output_path, kmer_size=20, threshold=60):
    query_sequences = parse_sequence_files(query_path)
    background_sequences = parse_sequence_files(background_path)

    output_records = annotate_homology(query_sequences, background_sequences, kmer_size, threshold)

    output_format = (
        "genbank" if output_path.endswith((".gb", ".gbk", ".genbank")) else "fasta"
    )
    with open(output_path, "w") as handle:
        SeqIO.write(output_records, handle, output_format)

def annotate_homology(query_sequences, background_sequences, kmer_size=20, threshold=60):
    bg_trie = process_background_sequences(background_sequences, kmer_size, merge = False)
    query_tries = process_query_sequences(query_sequences, kmer_size)
    output_records = []
    for query_seq in query_sequences:
        regions = find_non_homologous_regions(
            query_seq, query_tries[query_seq.id], bg_trie, kmer_size, threshold
        )
        output_record = create_annotated_record(query_seq, regions)
        output_records.append(output_record)
    
    return output_records

    