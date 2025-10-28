# %%
import pickle
import gc
import marisa_trie as mt
import tracemalloc
import time
import psutil

from choppy.homology_finder import parse_sequence_files, process_background_sequences

# %%
sequences = parse_sequence_files("data/S_cerevisiae-R64-GCA_000146045_cat.fa")
trie = process_background_sequences(sequences, kmer_size=25)
# %%
with open("data/saved_trie_25.pkl", "wb") as f:
    pickle.dump(trie, f)

# %%
keys = set(trie.keys())
print(type(keys))
with open("data/set_keys.pkl", "wb") as f:
    pickle.dump(keys, f)

# %%
del trie
gc.collect()

tracemalloc.start()
trie = process_background_sequences(sequences, kmer_size=20)
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
print(f"Peak memory: {peak / 1024 / 1024:.2f} MB")

# %%
trie.save("data/marisa_trie_20.marisa")

# %%
loaded_trie = mt.Trie()

del loaded_trie
loaded_trie = mt.Trie()
gc.collect()

tracemalloc.start()
loaded_trie.load("data/marisa_trie_20.marisa")
time.sleep(10)
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
print(f"Peak memory: {peak / 1024 / 1024:.2f} MB")

# %%
print(f"Loaded trie has {len(loaded_trie)} keys.")

# %%
loaded_trie = mt.Trie()

del loaded_trie
loaded_trie = mt.Trie()
gc.collect()
time.sleep(5)

process = psutil.Process()
mem_before = process.memory_info().rss
loaded_trie.load("data/marisa_trie_20.marisa")
time.sleep(10)
mem_after = process.memory_info().rss
print("Memory increase:", (mem_after - mem_before) / 1024 / 1024, "MB")

# %%
keys = set()

del keys
gc.collect()
time.sleep(5)

process = psutil.Process()
mem_before = process.memory_info().rss
keys = set(loaded_trie.keys())
time.sleep(10)
mem_after = process.memory_info().rss
print("Memory increase:", (mem_after - mem_before) / 1024 / 1024, "MB")
