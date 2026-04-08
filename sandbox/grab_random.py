# %%
import random
import pandas as pd
from Bio import SeqIO

# %%
# Let's grab sequences around genes since other regions mignt be strange
fasta_path = "/home/tyranchick/sds/sd17l002/p/dcm_mouse/data/reference/Mus_musculus.GRCm39.dna.primary_assembly.fa"
output_path = "data/random_seqs_50k.fa"
genome = {record.id: record.seq for record in SeqIO.parse(fasta_path, "fasta")}

# %%
gene_list = pd.read_csv("/home/tyranchick/sds/sd17l002/p/dcm_mouse/new_dcm_vault/data/reference/all_genes_GRCm39.107.txt")
gene_list['chr'] = gene_list['chr'].astype(str)
gene_list = gene_list[gene_list['chr'].isin([str(i) for i in range(1, 20)] + ['X', 'Y'])]

# %%
sequences = {}
n = 20
reg_len = 50000

for seq_id in range(n):
    gene = gene_list.sample(1).iloc[0]
    random_pos = random.randint(gene['start'] - reg_len, gene['end'])
    seq = genome.get(gene['chr'], "")[random_pos:random_pos + reg_len]
    sequences[f"seq_{seq_id}"] = seq

with open(output_path, 'w') as f:
    for idx in sequences:
        f.write(f">{idx}\n{sequences[idx]}\n")

