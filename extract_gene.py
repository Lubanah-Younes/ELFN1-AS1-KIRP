# extract_genes.py
import pandas as pd
import numpy as np
from scipy import stats
import gzip
import os

print("Extracting top 200 genes correlated with ELFN1-AS1...")

# Load data (assuming files are in current directory)
GENE = "ENSG00000238117"
EXPR_FILE = "tcga_RSEM_gene_tpm.gz"
CLIN_FILE = "Survival_SupplementalTable_S1_20171025_xena_sp"

# Load KIRP samples
clin = pd.read_csv(CLIN_FILE, sep='\t')
kirp_samples = clin[clin['cancer type abbreviation'] == 'KIRP']['sample'].tolist()
print(f"Found {len(kirp_samples)} KIRP samples")

# Load expression data
with gzip.open(EXPR_FILE, 'rt') as f:
    header = f.readline().strip().split('\t')
    kirp_indices = [i for i, s in enumerate(header[1:]) if s in kirp_samples]
    kirp_sample_names = [header[i+1] for i in kirp_indices]
    
    # Find gene
    gene_found = False
    for line in f:
        parts = line.strip().split('\t')
        if parts[0].startswith(GENE):
            gene_expr = [float(parts[i+1]) for i in kirp_indices]
            gene_found = True
            break
    
    if not gene_found:
        print("Gene not found!")
        exit()

# Calculate correlations with all genes (simplified)
print("Calculating correlations...")
correlations = []
gene_idx = 0

# We need to read the file again
with gzip.open(EXPR_FILE, 'rt') as f:
    f.readline()  # skip header
    for line in f:
        parts = line.strip().split('\t')
        if parts[0].startswith(GENE):
            continue  # skip our gene
        
        if gene_idx % 1000 == 0:
            print(f"Processed {gene_idx} genes...")
        
        other_expr = [float(parts[i+1]) for i in kirp_indices]
        corr, p = stats.pearsonr(gene_expr, other_expr)
        correlations.append({'gene': parts[0], 'correlation': corr, 'p_value': p})
        gene_idx += 1

# Create DataFrame and save
corr_df = pd.DataFrame(correlations).sort_values('correlation', ascending=False)
corr_df.to_csv("ELFN1-AS1_correlations.csv", index=False)

# Save top 200 positive
top_pos = corr_df.head(200)['gene'].tolist()
with open("ELFN1-AS1_top200_positive_genes.txt", "w") as f:
    f.write("\n".join(top_pos))
print("✅ Top 200 positive genes saved to ELFN1-AS1_top200_positive_genes.txt")

# Save top 200 negative
top_neg = corr_df.tail(200)['gene'].tolist()
with open("EFLFN1-AS1_top200_negative_genes.txt", "w") as f:
    f.write("\n".join(top_neg))
print("✅ Top 200 negative genes saved to ELFN1-AS1_top200_negative_genes.txt")