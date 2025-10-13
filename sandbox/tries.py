# %%
import marisa_trie as mt
import sys
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