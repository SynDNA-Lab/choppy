# %%
from Bio import SeqIO
from choppy.homology_finder import load_trie


# %%
sequences = list(SeqIO.parse("data/random_seqs_50k.fa", "fasta"))

test_seq = str(sequences[0].seq)
kmer_size = 15
bg_trie = load_trie("../data/S_cerevisiae-R64-GCA_000146045_cat_15.marisa")


