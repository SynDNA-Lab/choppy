# %%
import marisa_trie as mt
import sys
import tracemalloc
import gc
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
sys.path.insert(0, '/home/tyranchick/git/TARnche')
from tarnche.homology_finder import (_parse_sequence_files, create_kmer_dictionary, create_kmer_trie)
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
# %%
type(trie)
# %%
trie = create_kmer_trie(seqs[0], 15, True)

print(set(trie.keys()) == set(kmer_dict.keys()))

# %%
extra_keys = set(trie.keys()) - set(kmer_dict.keys())
print("Keys in trie but not in kmer_dict:", extra_keys)
print(str(seqs[0].seq).find(list(extra_keys)[1]))
seqs[0].seq[48487:48487+15]