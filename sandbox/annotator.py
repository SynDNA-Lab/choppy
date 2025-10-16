# %%
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from tarnche.fragment_annotator import extract_no_homology_regions
# %%

input_record = list(SeqIO.parse("/home/tyranchick/git/TARnche/tests/data/annotated.gb", "genbank"))
nohom_regions = extract_no_homology_regions(input_record[1])
nohom_regions
# %%
input_record[1].features[0]["label"]