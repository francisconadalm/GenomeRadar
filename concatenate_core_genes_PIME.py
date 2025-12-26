import pandas as pd
from collections import defaultdict, OrderedDict

df = pd.read_csv("table_proteins_PIME.txt", sep="\t")

gene_order = ['Integrase', 'AlpA', 'Rep', 'TrbJ', 'TrbL', 'LuxR', 'TraJ', 'TraI']

grouped = defaultdict(dict)
for _, row in df.iterrows():
    region = row['seq_name']
    gene_type = row['gene_name']
    translation = row['Translation']
    if gene_type not in grouped[region]:
        grouped[region][gene_type] = translation

with open("concatented_PIME_core.fasta", "w") as fasta:
    for region, genes in grouped.items():
        sequence_parts = [genes[gene] for gene in gene_order if gene in genes]
        if sequence_parts:
            sequence = ''.join(sequence_parts)
            fasta.write(f">{region}\n")
            for i in range(0, len(sequence), 60):
                fasta.write(sequence[i:i+60] + "\n")
