# %%
import marisa_trie as mt
import sys
import tracemalloc
import gc
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
sys.path.insert(0, '/home/tyranchick/git/TARnche')
from tarnche.homology_finder import _parse_sequence_files
import pickle
# %%
keys = ['foo', 'bar', 'foobar', 'foo']
values = [(1, 2), (2, 1), (3, 3), (2, 1)]
fmt = "<HH"   # two short integers.
trie = mt.RecordTrie(fmt, zip(keys, values))
sys.getsizeof(trie)
# %%
d = dict(zip(keys, values))
sys.getsizeof(d)
# %%
trie['foo']
# %%
b_trie = mt.RecordTrie('<B', [('AAA', (1,)), ('CCC', (0,)), ('AAA', (0,))])
b_trie['AAA']
# %%
seqs = _parse_sequence_files('data/lambda_genome-NC_001416.gb')

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

kmer_dict = create_kmer_dictionary(seqs[0], 15)
sys.getsizeof(kmer_dict) / 1024 / 1024
# %%
len(kmer_dict.keys())
# %%
del trie
gc.collect()

tracemalloc.start()
trie = mt.Trie(kmer_dict.keys())
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
print(f"Peak memory: {peak / 1024 / 1024:.2f} MB")