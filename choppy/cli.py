import click

@click.group()
def cli():
    """Choppy: Tools for chopping DNA sequences for TAR cloning."""
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
    from choppy.homology_finder import file_process_homology
    
    click.echo(f"Processing query sequence: {query_sequence}")
    click.echo(f"Background sequences: {', '.join(background)}")
    click.echo(f"K-mer size: {kmer_size}, Threshold: {threshold}")
    
    # Call the main function
    file_process_homology(
        query_path=query_sequence,
        background_path=list(background),
        output_path=output,
        kmer_size=kmer_size,
        threshold=threshold
    )
    
    click.echo(f"✓ Wrote annotated sequences to {output}")


@cli.command(
    name="fragment",
    help="""
    Fragment an annotated sequence around no-homology regions with overlap and boundary constraints.
    INPUT is the GenBank file with no-homology annotated regions (from annotate-homology command).
    """,
    short_help="Fragment a sequence with overlap constraints.",
    no_args_is_help=True,
)
@click.argument(
    "input",
    type=click.Path(exists=True, readable=True),
)
@click.option(
    "--min-size",
    type=int,
    required=True,
    help="Minimum fragment size (bp)",
)
@click.option(
    "--max-size",
    type=int,
    required=True,
    help="Maximum fragment size (bp)",
)
@click.option(
    "--min-overlap",
    type=int,
    required=True,
    help="Minimum overlap between fragments (bp)",
)
@click.option(
    "--max-overlap",
    type=int,
    required=True,
    help="Maximum overlap between fragments (bp)",
)
@click.option(
    "--boundary-motif",
    type=str,
    default=None,
    help="Motif that must occur at fragment boundaries (start and end)",
)
@click.option(
    "--min-step",
    type=int,
    default=10,
    show_default=True,
    help="Step size for scanning overlap positions (bp)",
)
@click.option(
    "-o",
    "--output",
    default="fragmented.gb",
    show_default=True,
    type=click.Path(writable=True),
    help="Output annotated GenBank file",
)
@click.option(
    "-f",
    "--fasta",
    default="fragments.fasta",
    show_default=True,
    type=click.Path(writable=True),
    help="Output multi-fasta file for fragments",
)
def fragment_cmd(input, min_size, max_size, min_overlap, max_overlap, boundary_motif, min_step, output, fasta):
    """Fragment a sequence with overlap constraints."""
    from choppy.fragment_annotator import fragment_from_file, FragmentConfig
    
    click.echo(f"Processing input file: {input}")
    click.echo(f"Fragment size: {min_size}-{max_size} bp")
    click.echo(f"Overlap: {min_overlap}-{max_overlap} bp")
    if boundary_motif:
        click.echo(f"Boundary motif: {boundary_motif}")
    click.echo(f"Step size: {min_step} bp")
    
    # Create configuration
    config = FragmentConfig(
        min_size=min_size,
        max_size=max_size,
        min_overlap=min_overlap,
        max_overlap=max_overlap,
        motif=boundary_motif,
        min_step=min_step,
    )
    
    # Call the fragmentation function
    fragment_from_file(
        input_file=input,
        output_gb=output,
        output_fasta=fasta,
        config=config,
    )
    
    click.echo(f"✓ Fragmentation complete!")


@cli.command(
    name="store-trie",
    help="""
    Create and store a k-mer trie from sequence file(s) for later use.
    INPUT is the path to sequence file(s) (GenBank or FASTA).
    The trie will be saved to the specified output file for use with other commands.
    """,
    short_help="Create and store k-mer trie from sequences.",
    no_args_is_help=True,
)
@click.argument(
    "input",
    type=click.Path(exists=True, readable=True),
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
    "-o",
    "--output",
    required=True,
    type=click.Path(writable=True),
    help="Output trie file path",
)
def store_trie_cmd(input, kmer_size, output):
    """Create and store a k-mer trie from sequence file(s)."""
    from choppy.homology_finder import store_trie_from_file
    
    click.echo(f"Processing input file: {input}")
    click.echo(f"K-mer size: {kmer_size}")
    
    # Call the store_trie_from_file function
    store_trie_from_file(
        path=input,
        trie_path=output,
        kmer_size=kmer_size,
    )
    
    click.echo(f"✓ Trie saved to {output}")


