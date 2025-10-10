# TARnche
Functionalities for chopping up DNA for TAR cloning (and other homology based methods).
## Motivation
This repo is for developing tools to chop DNA sequences into smaller fragments
for in vivo DNA assembly by homologous recombination in the yeast S. cerevisiae
(TAR cloning).

**Name:** the name is inspired from the term 'la tranche' (frz. for slice) and mixed
with the term 'TAR' (tranformation associated recombination**.

## Usage
### 0) Set-Up
 - requirements (Python packages):
  - Bio
  - tqdm

### 1) Annotate homology free regions
```bash
python3 homology_finder.py lambda_genome-NC_001416.gb \
    -b S_cerevisiae-R64-GCA_000146045_cat.fa pSDL32-TAR_shuffle-pCC1-Gibson.gb \
    pMaM819-Double-marker-plasmid_fixedColEI-pMLPstar.gb \
    -k 15 -t 60 -o lambda_genome-NC_001416_hfanno.gb
```
### 2) Fragment the DNA sequence based on homology free regions and further constraints
```bash
python3 fragment_annotator.py lambda_genome-NC_001416_hfanno.gb \
    --min-size 500 --max-size 3188 --min-overlap 60 --max-overlap 100  \
    --boundary-motif 'GC' -o lambda_genome-NC_001416_hfanno_fragments.gb \
    -f lambda_genome-NC_001416_hfanno_fragments.fa
```
## Acknowledgement
Initial prototype (not part of this code base) was written by [@dglaizal](https://github.com/dglaizal)
