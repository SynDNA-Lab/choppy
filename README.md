# TARnche
Functionalities for chopping up DNA for TAR cloning (and other homology based methods).

## Motivation
This repo is for developing tools to chop DNA sequences into smaller fragments
for in vivo DNA assembly by homologous recombination in the yeast S. cerevisiae
(TAR cloning).

**Name:** the name is inspired from the term 'la tranche' (frz. for slice) and mixed
with the term 'TAR' (transformation associated recombination).

## Installation

### From GitHub repository

```bash
# Clone the repository
git clone https://github.com/SynDNA-Lab/TARnche.git
cd TARnche

# Install with poetry (recommended)
poetry install

# Or install with pip
pip install .
```

## Usage

### 1) Annotate homology-free regions

Find and annotate non-homologous regions in query sequences compared to background sequences:

```bash
tarnche annotate-homology lambda_genome-NC_001416.gb \
    -b S_cerevisiae-R64-GCA_000146045_cat.fa \
    -b pSDL32-TAR_shuffle-pCC1-Gibson.gb \
    -b pMaM819-Double-marker-plasmid_fixedColEI-pMLPstar.gb \
    -k 15 -t 60 -o lambda_genome-NC_001416_hfanno.gb
```

**Options:**
- `QUERY_SEQUENCE` - Path to the query sequence file (GenBank or FASTA)
- `-b, --background` - Path(s) to background sequence file(s) (can be specified multiple times)
- `-k, --kmer-size` - Size of k-mers (default: 20)
- `-t, --threshold` - Minimum size of non-homologous regions to report (default: 60)
- `-o, --output` - Output file name (default: annotated_sequences.gb)

### 2) Fragment the DNA sequence

Fragment an annotated sequence based on homology-free regions and overlap constraints:

```bash
tarnche fragment lambda_genome-NC_001416_hfanno.gb \
    --min-size 500 --max-size 3188 \
    --min-overlap 60 --max-overlap 100 \
    --boundary-motif 'GC' \
    -o lambda_genome-NC_001416_hfanno_fragments.gb \
    -f lambda_genome-NC_001416_hfanno_fragments.fa
```

**Options:**
- `INPUT` - GenBank file with annotated no-homology regions (from step 1)
- `--min-size` - Minimum fragment size in bp (required)
- `--max-size` - Maximum fragment size in bp (required)
- `--min-overlap` - Minimum overlap between fragments in bp (required)
- `--max-overlap` - Maximum overlap between fragments in bp (required)
- `--boundary-motif` - Motif that must occur at fragment boundaries (optional)
- `--min-step` - Step size for scanning overlap positions in bp (default: 10)
- `-o, --output` - Output annotated GenBank file (default: fragmented.gb)
- `-f, --fasta` - Output multi-fasta file for fragments (default: fragments.fasta)

### Help

For more information on any command:

```bash
tarnche --help
tarnche annotate-homology --help
tarnche fragment --help
```

## Acknowledgement
Initial prototype (not part of this code base) was written by [@dglaizal](https://github.com/dglaizal)

