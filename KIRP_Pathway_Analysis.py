# =============================================================================
# Pathway Analysis for ELFN1-AS1 in KIRP
# Author: Lubanah Younes
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gseapy as gp
from scipy import stats
import gzip
import os

print("="*50)
print("KIRP PATHWAY ANALYSIS FOR ELFN1-AS1")
print("="*50)

# =============================================================================
# 1. Load KIRP sample list
# =============================================================================
print("\n[1/5] Loading KIRP sample list...")
clin = pd.read_csv("Survival_SupplementalTable_S1_20171025_xena_sp", sep='\t')
kirp_samples = clin[clin['cancer type abbreviation'] == 'KIRP']['sample'].tolist()
print(f"✅ Found {len(kirp_samples)} KIRP samples")

# =============================================================================
# 2. Locate ELFN1-AS1 and extract its expression for KIRP samples
# =============================================================================
print("\n[2/5] Locating ELFN1-AS1...")
gene_target = "ENSG00000238117"
gene_id = None
gene_expr = None
kirp_indices = []
kirp_sample_names = []

with gzip.open("tcga_RSEM_gene_tpm.gz", 'rt') as f:
    header = f.readline().strip().split('\t')
    all_samples = header[1:]
    
    # Identify KIRP sample indices and names
    kirp_indices = [i for i, s in enumerate(all_samples) if s in kirp_samples]
    kirp_sample_names = [all_samples[i] for i in kirp_indices]
    print(f"✅ Found {len(kirp_sample_names)} KIRP samples in expression file")
    
    # Find the target gene
    found_gene = False
    for line_num, line in enumerate(f):
        parts = line.strip().split('\t')
        current_gene = parts[0]
        
        if current_gene.startswith(gene_target):
            gene_id = current_gene
            # Extract values for KIRP samples
            values = [float(parts[i+1]) for i in kirp_indices]
            gene_expr = pd.Series(values, index=kirp_sample_names)
            print(f"✅ Found gene at line {line_num}: {gene_id}")
            found_gene = True
            break
    
    if not found_gene:
        print("❌ Target gene not found")
        exit()

# =============================================================================
# 3. Calculate correlations with all other genes
# =============================================================================
print("\n[3/5] Calculating correlations (this takes 5-10 minutes)...")

correlations = []
total_genes = 0

# Reset file pointer (need to reopen)
with gzip.open("tcga_RSEM_gene_tpm.gz", 'rt') as f:
    # Skip header
    f.readline()
    
    for line_num, line in enumerate(f):
        if line_num % 5000 == 0:
            print(f"   Processed {line_num} genes...")
        
        parts = line.strip().split('\t')
        current_gene = parts[0]
        
        # Skip if this is our target gene
        if current_gene == gene_id:
            continue
        
        # Extract expression for KIRP samples
        values = [float(parts[i+1]) for i in kirp_indices]
        other_expr = pd.Series(values, index=kirp_sample_names)
        
        # Calculate correlation
        corr, p_val = stats.pearsonr(gene_expr, other_expr)
        correlations.append({
            'gene': current_gene,
            'correlation': corr,
            'p_value': p_val
        })
        total_genes += 1

print(f"✅ Calculated correlations for {total_genes} genes")

# Create DataFrame
corr_df = pd.DataFrame(correlations).sort_values('correlation', ascending=False)
print(f"✅ Correlation matrix shape: {corr_df.shape}")

# =============================================================================
# 4. Extract top correlated genes
# =============================================================================
print("\n[4/5] Extracting top correlated genes...")

top_pos = corr_df.head(200)['gene'].tolist()
top_neg = corr_df.tail(200)['gene'].tolist()
print("\n🔥 Top 10 positively correlated genes:")
print(corr_df.head(10)[['gene', 'correlation']].to_string())

print("\n🔥 Top 10 negatively correlated genes:")
print(corr_df.tail(10)[['gene', 'correlation']].to_string())

# =============================================================================
# 5. Pathway analysis (corrected)
# =============================================================================
print("\n[5/5] Running pathway enrichment analysis...")
output_dir = "results_KIRP_ENSG00000238117"
os.makedirs(output_dir, exist_ok=True)

# Remove NaN values from top correlated lists
top_pos_clean = [g for g in top_pos[:100] if not pd.isna(g)]
top_neg_clean = [g for g in top_neg[:100] if not pd.isna(g)]

print(f"   Positive genes after cleaning: {len(top_pos_clean)}")
print(f"   Negative genes after cleaning: {len(top_neg_clean)}")

try:
    # Positive correlations
    print("\n🔬 Pathways enriched with ELFN1-AS1 (positively correlated):")
    enr_pos = gp.enrichr(gene_list=top_pos_clean, 
                         gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2021'],
                         organism='human')  # Changed from 'Human' to 'human'
    print(enr_pos.results.head(10)[['Term', 'P-value']].to_string())
    enr_pos.results.to_csv(f"{output_dir}/KIRP_positive_pathways.csv", index=False)
    
    # Negative correlations
    print("\n🔬 Pathways suppressed by ELFN1-AS1 (negatively correlated):")
    enr_neg = gp.enrichr(gene_list=top_neg_clean, 
                         gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2021'],
                         organism='human')  # Changed from 'Human' to 'human'
    print(enr_neg.results.head(10)[['Term', 'P-value']].to_string())
    enr_neg.results.to_csv(f"{output_dir}/KIRP_negative_pathways.csv", index=False)
    
    print(f"\n✅ Results saved to {output_dir}")
except Exception as e:
    print(f"Pathway analysis error: {e}")

print("\n" + "="*50)
print("✅ ANALYSIS COMPLETE")
print("="*50)