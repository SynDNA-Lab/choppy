#!/usr/bin/env python
# coding: utf-8

import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from tqdm import tqdm
from typing import List, Dict, Union


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
        format_type = "genbank" if file_ext in [".gb", ".gbk", ".genbank"] else "fasta"

        if file_ext not in [".gb", ".gbk", ".genbank", ".fa"]:
            raise Exception(f"Filetype of file '{file_path}' is not supported.")

        # Parse all sequences in the file
        records.extend(list(SeqIO.parse(file_path, format_type)))

    return records


def create_kmer_dictionary(sequence: SeqRecord, kmer_size: int) -> dict:
    """
    Create a dictionary of k-mers and their positions in the sequence.

    Args:
        sequence (SeqRecord): DNA sequence record
        kmer_size (int): Size of k-mers

    Returns:
        dict: Dictionary with k-mers as keys and positions as values
    """
    kmer_dict = {}

    # Check if sequence is circular
    is_circular = sequence.annotations.get("topology", "").lower() == "circular"

    # Prepare sequence for processing
    sequence_str = str(sequence.seq)
    if is_circular:
        sequence_larger = sequence_str + sequence_str[0:kmer_size]
    else:
        sequence_larger = sequence_str

    # Forward strand
    for nucleotide in tqdm(
        range(1, len(sequence_str) + 1), desc="Processing forward strand"
    ):
        if nucleotide + kmer_size <= len(sequence_larger):
            kmer = sequence_larger[nucleotide - 1 : nucleotide - 1 + kmer_size]
            if kmer not in kmer_dict:
                kmer_dict[kmer] = [nucleotide]
            else:
                kmer_dict[kmer].append(nucleotide)

    # Reverse complement strand
    sequence_rc = sequence.seq.reverse_complement()
    sequence_rc_str = str(sequence_rc)
    if is_circular:
        sequence_rc_larger = sequence_rc_str + sequence_rc_str[0:kmer_size]
    else:
        sequence_rc_larger = sequence_rc_str

    for nucleotide in tqdm(
        range(1, len(sequence_rc_str) + 1), desc="Processing reverse strand"
    ):
        if nucleotide + kmer_size <= len(sequence_rc_larger):
            kmer = sequence_rc_larger[nucleotide - 1 : nucleotide - 1 + kmer_size]
            if kmer not in kmer_dict:
                kmer_dict[kmer] = [-nucleotide]
            else:
                kmer_dict[kmer].append(-nucleotide)

    return kmer_dict


def update_kmer_dictionary(
    kmer_dict: dict, sequence: SeqRecord, kmer_size: int, offset: int
):
    """
    Update existing k-mer dictionary with positions from another sequence.

    Args:
        kmer_dict (dict): Existing k-mer dictionary
        sequence (SeqRecord): DNA sequence record to process
        kmer_size (int): Size of k-mers
        offset (int): Position offset for the new sequence
    """
    # Check if sequence is circular
    is_circular = sequence.annotations.get("topology", "").lower() == "circular"

    # Forward strand
    sequence_str = str(sequence.seq)
    if is_circular:
        sequence_larger = sequence_str + sequence_str[0:kmer_size]
    else:
        sequence_larger = sequence_str

    for nucleotide in tqdm(
        range(1, len(sequence_str) + 1), desc="Updating forward strand"
    ):
        if nucleotide + kmer_size <= len(sequence_larger):
            kmer = sequence_larger[nucleotide - 1 : nucleotide - 1 + kmer_size]
            if kmer in kmer_dict:
                kmer_dict[kmer].append(nucleotide + offset)

    # Reverse complement strand
    sequence_rc = sequence.seq.reverse_complement()
    sequence_rc_str = str(sequence_rc)
    if is_circular:
        sequence_rc_larger = sequence_rc_str + sequence_rc_str[0:kmer_size]
    else:
        sequence_rc_larger = sequence_rc_str

    for nucleotide in tqdm(
        range(1, len(sequence_rc_str) + 1), desc="Updating reverse strand"
    ):
        if nucleotide + kmer_size <= len(sequence_rc_larger):
            kmer = sequence_rc_larger[nucleotide - 1 : nucleotide - 1 + kmer_size]
            if kmer in kmer_dict:
                kmer_dict[kmer].append(-nucleotide - offset)


def process_sequences(
    background_sequences: List[SeqRecord],
    query_sequences: List[SeqRecord],
    kmer_size: int,
) -> Dict[str, dict]:
    """
    Process background and query sequences to create k-mer dictionaries.

    Args:
        background_sequences (List[SeqRecord]): List of background sequences
        query_sequences (List[SeqRecord]): List of query sequences
        kmer_size (int): Size of k-mers

    Returns:
        Dict[str, dict]: Dictionary of k-mer dictionaries for each query sequence
    """
    results = {}

    for query_seq in query_sequences:
        # Initialize empty k-mer dictionary
        kmer_dict = {}

        # Process query sequences first
        offset = 0
        kmer_dict = create_kmer_dictionary(query_seq, kmer_size)
        offset += len(query_seq)
        # Update dictionary with background sequences
        for bg_seq in background_sequences:
            update_kmer_dictionary(kmer_dict, bg_seq, kmer_size, offset)
            offset += len(bg_seq)

        # Store results for this query sequence
        results[query_seq.id] = kmer_dict

    return results


def find_non_homologous_regions(
    sequence_length: int, kmer_dict: dict, kmer_size: int, threshold: int
) -> dict:
    """
    Find regions without homology based on k-mer positions.

    Args:
        sequence_length (int): Length of the sequence
        kmer_dict (dict): Dictionary of k-mers and their positions
        kmer_size (int): Size of k-mers
        threshold (int): Minimum size of non-homologous regions to report

    Returns:
        dict: Dictionary of non-homologous regions
    """
    # Create a list of all nucleotide positions
    nucleotide_positions = [i for i in range(sequence_length)]

    # Mark homologous positions
    for kmer, positions in kmer_dict.items():
        if len(positions) > 1:
            for pos in positions:
                if sequence_length >= pos > 0:
                    nucleotide_positions[pos : pos + kmer_size] = ["h"] * kmer_size
                elif -sequence_length <= pos < 0:
                    start = pos + 1 + sequence_length - kmer_size
                    end = pos + 1 + sequence_length
                    nucleotide_positions[start:end] = ["h"] * kmer_size

    # Filter out homologous positions
    non_homologous = [i for i in nucleotide_positions if i != "h"]

    # Find continuous non-homologous regions
    regions = {}
    region_count = 1
    i = 0

    while i < len(non_homologous):
        start = non_homologous[i]
        while (
            i + 1 < len(non_homologous)
            and non_homologous[i] + 1 == non_homologous[i + 1]
        ):
            i += 1
        end = non_homologous[i]

        if end - start >= threshold:
            regions[f"Region_{region_count}"] = (start, end)
            region_count += 1
        i += 1

    return regions


def create_annotated_record(
    sequence: SeqRecord, regions: dict, output_format: str = "genbank"
) -> SeqRecord:
    """
    Create an annotated sequence record with non-homologous regions.

    Args:
        sequence (SeqRecord): Input sequence record
        regions (dict): Dictionary of non-homologous regions
        output_format (str): Output file format ('genbank' or 'fasta')

    Returns:
        SeqRecord: Annotated sequence record
    """
    features = getattr(sequence, "features", [])

    # Add features for non-homologous regions
    for region_name, (start, end) in regions.items():
        feature = SeqFeature(
            FeatureLocation(start, end),
            type="misc_feature",
            qualifiers={
                "note": f"Non-homologous region {region_name}",
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


def main():
    parser = argparse.ArgumentParser(
        description="Find non-homologous regions between sequences"
    )
    parser.add_argument(
        "query_sequence", help="Path to query sequence file(s) (GenBank or FASTA)"
    )
    parser.add_argument(
        "-b",
        "--background",
        nargs="+",
        required=True,
        help="Path(s) to background sequence file(s) (GenBank or FASTA)",
    )
    parser.add_argument(
        "-k", "--kmer-size", type=int, default=20, help="Size of k-mers (default: 20)"
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=int,
        default=60,
        help="Minimum size of non-homologous regions to report (default: 60)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="annotated_sequences.gb",
        help="Output file name (default: annotated_sequences.gb)",
    )

    args = parser.parse_args()

    # Parse input sequences
    query_sequences = parse_sequence_files(args.query_sequence)
    background_sequences = parse_sequence_files(args.background)

    # Process sequences and get k-mer dictionaries
    kmer_dicts = process_sequences(
        background_sequences, query_sequences, args.kmer_size
    )

    # Process each query sequence
    output_records = []
    for query_seq in query_sequences:
        # Find non-homologous regions
        regions = find_non_homologous_regions(
            len(query_seq.seq), kmer_dicts[query_seq.id], args.kmer_size, args.threshold
        )

        # Create annotated record
        output_record = create_annotated_record(query_seq, regions)
        output_records.append(output_record)

    # Write output file
    output_format = (
        "genbank" if args.output.endswith((".gb", ".gbk", ".genbank")) else "fasta"
    )
    with open(args.output, "w") as handle:
        SeqIO.write(output_records, handle, output_format)

    print(
        f"Processed {len(query_sequences)} query sequences against {len(background_sequences)} background sequences"
    )
    print(f"Wrote results to {args.output}")


if __name__ == "__main__":
    main()
