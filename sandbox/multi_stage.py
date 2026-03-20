# %%
from Bio import SeqIO
from collections import Counter

# %%
sequences = list(SeqIO.parse("data/random_seqs_50k.fa", "fasta"))
test_seq = str(sequences[0].seq)

test_seq = "ATGCATGCATGCATGC"

neigh_size = 8
kmer_size = 4

kmers = dict()
total_kmers = len(test_seq) - kmer_size + 1
hom_free = [0] * total_kmers

i = 0
j = 0
while j < total_kmers:
    if i < total_kmers:
        kmer = test_seq[i:i + kmer_size]
        last_pos = kmers.get(kmer)
        if last_pos is None or i + kmer_size - last_pos > neigh_size:
            print("no previous occurrence at position", i)
            hom_free[i] += 1
        
        kmers[kmer] = i

    if i >= neigh_size - kmer_size:
        kmer = test_seq[j:j + kmer_size]
        last_pos = kmers.get(kmer)
        if last_pos == j:
            print("no future occurrence at position", j)
            hom_free[j] += 1
        j += 1
    print(kmers, hom_free)
    i += 1

#hom_free.extend([2] * (kmer_size - 1))

# %%


threshold = 100
regions = []
start = None

for i, val in enumerate(hom_free):
    if val == 2:
        if start is None:
            start = i
    else:
        if start is not None and i - start + kmer_size - 1 >= threshold:
            regions.append((start, i + kmer_size - 1))
        start = None

if start is not None and len(hom_free) - start >= threshold:
    regions.append((start, len(hom_free) + kmer_size - 1))

