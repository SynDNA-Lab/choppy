import click
from Bio import SeqIO


@click.group()
def cli():
    """TARnche: Tools for chopping DNA sequences for TAR cloning."""
    pass


@cli.command(
    name="annotate-homology",
    help="""
    Find non-homologous regions in query sequences compared to background sequences and itself.
    QUERY_SEQUENCE is the path to the query sequence file (GenBank or FASTA).
    Outputs annotated GenBank files with non-homologous regions marked as features.
    """,
    short_help="Find and annotate non-homologous regions in sequences.",
    no_args_is_help=True,

)
@click.argument(
    "query_sequence",
    type=click.Path(exists=True, readable=True),
)
@click.option(
    "-b",
    "--background",
    multiple=True,
    type=click.Path(exists=True, readable=True),
    help="Path(s) to background sequence file(s) (GenBank or FASTA)",
)
@click.option(
    "-k",
    "--kmer-size",
    type=int,
    default=20,
    show_default=True,
    help="Size of k-mers",
)
@click.option(
    "-t",
    "--threshold",
    type=int,
    default=60,
    show_default=True,
    help="Minimum size of non-homologous regions to report",
)
@click.option(
    "-o",
    "--output",
    default="annotated_sequences.gb",
    show_default=True,
    type=click.Path(writable=True),
    help="Output file name",
)
def annotate_homology_cmd(query_sequence, background, kmer_size, threshold, output):
    """Find and annotate non-homologous regions."""
    from tarnche.homology_finder import annotate_homology
    
    click.echo(f"Processing query sequence: {query_sequence}")
    click.echo(f"Background sequences: {', '.join(background)}")
    click.echo(f"K-mer size: {kmer_size}, Threshold: {threshold}")
    
    # Call the main function
    annotate_homology(
        query_path=query_sequence,
        background_path=list(background),
        output_path=output,
        kmer_size=kmer_size,
        threshold=threshold
    )
    
    click.echo(f"✓ Wrote annotated sequences to {output}")
