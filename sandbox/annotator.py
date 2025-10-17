# %%
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import networkx as nx

from tarnche.fragment_annotator import extract_no_homology_regions, Fragmentor, FragmentConfig
# %%

input_record = list(SeqIO.parse("/home/tyranchick/git/TARnche/tests/data/annotated.gb", "genbank"))
nohom_regions = extract_no_homology_regions(input_record[1])
nohom_regions
# %%
input_record[1].features[0]["label"]
# %%
seq = "A" * 1000
nohom_regions = [(50, 100), (800, 850)]
config = FragmentConfig(100, 200, 20, 50, min_step=30)
        
fragmentor = Fragmentor(seq, nohom_regions, config)
fragmentor.graph.edges
# %%
fragmentor.fix_graph()
fragmentor.graph
# %%
seq = "A" * 100 + "GAATTC" + "A" * 20 + "GAATTC" + "A" * 400 + "GAATTC" + "A" * 20 + "GAATTC" + "A" * 100
nohom_regions = [(90, 150), (200, 400), (500, len(seq))]
config = FragmentConfig(100, 250, 20, 50, min_step=20, motif="GAATTC")
        
fragmentor = Fragmentor(seq, nohom_regions, config)
initial_nodes = list(fragmentor.graph.nodes)
initial_nodes
# %%        
components = list(nx.weakly_connected_components(fragmentor.graph))
fragmentor.keep_homology_constraint(components)

fragmentor.graph.nodes

# %%
seq = "A" * 1000
nohom_regions = [(50, 100), (900, 950)]
config = FragmentConfig(100, 250, 20, 50, min_step=20)

fragmentor = Fragmentor(seq, nohom_regions, config)
initial_nodes = list(fragmentor.graph.nodes)
initial_nodes

# %%
seq = "A" * 1000
nohom_regions = [(50, 100), (800, 850)]
config = FragmentConfig(100, 450, 20, 50, min_step=10)

fragmentor = Fragmentor(seq, nohom_regions, config)
initial_nodes = list(fragmentor.graph.nodes)
print(initial_nodes)
print(list(fragmentor.graph.edges))
# %%

# Graph might be disconnected initially
path = fragmentor.get_shortest_path()
print(path)

# %%
seq = "A" * 650
nohom_regions = []
config = FragmentConfig(100, 350, 20, 50, min_step=10)

fragmentor = Fragmentor(seq, nohom_regions, config)
fragmentor.graph.nodes

# %%
seq = "A" * 300
nohom_regions = []
config = FragmentConfig(100, 250, 20, 50)

fragmentor = Fragmentor(seq, nohom_regions, config)
print(fragmentor.graph.nodes)
components = list(nx.weakly_connected_components(fragmentor.graph))
print(components)

fragmentor.fix_graph()
print(fragmentor.graph.nodes)

path = fragmentor.get_shortest_path()
print(path)