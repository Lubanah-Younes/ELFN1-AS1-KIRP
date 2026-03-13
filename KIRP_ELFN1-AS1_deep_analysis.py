import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
import gzip
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from lifelines import CoxPHFitter
import shap
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import os
import warnings
warnings.filterwarnings('ignore')

GENE = "ENSG00000238117"  # ELFN1-AS1
CANCER = "KIRP"
EXPR_FILE = "tcga_RSEM_gene_tpm.gz"
CLIN_FILE = "Survival_SupplementalTable_S1_20171025_xena_sp"
OUTPUT_DIR = f"results_{CANCER}_{GENE}"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("="*70)
print(f"🔬 DEEP ANALYSIS OF {GENE} IN {CANCER}")
print("="*70)

# 1. Load expression data
print("\n[1/5] Loading expression data...")
with gzip.open(EXPR_FILE, 'rt') as f:
    header = f.readline().strip().split('\t')
    for line in f:
        if line.startswith(GENE):
            parts = line.strip().split('\t')
            gene_id_full = parts[0]
            gene_expr_all = pd.Series(index=header[1:], data=parts[1:], dtype=float)
            break
    else:
        raise ValueError(f"Gene {GENE} not found")
print(f"✅ Found gene: {gene_id_full}")

# 2. Load clinical data
print("\n[2/5] Loading clinical data...")
clin = pd.read_csv(CLIN_FILE, sep='\t')
print(f"✅ Clinical columns: {clin.columns.tolist()}")

# 3. Subset KIRP samples
print("\n[3/5] Subsetting KIRP samples...")
clin_kirp = clin[clin['cancer type abbreviation'] == CANCER].copy()
print(f"✅ KIRP samples in clinical data: {clin_kirp.shape[0]}")

# 4. Align expression and clinical samples
common_samples = [s for s in clin_kirp['sample'] if s in gene_expr_all.index]
print(f"✅ Samples with both expression and clinical data: {len(common_samples)}")

expr_kirp = gene_expr_all[common_samples]
clin_kirp = clin_kirp[clin_kirp['sample'].isin(common_samples)]

# 5. Prepare survival data
print("\n[4/5] Preparing survival data...")
clin_clean = clin_kirp[['sample', 'OS', 'OS.time']].dropna().copy()
clin_clean['OS.time'] = pd.to_numeric(clin_clean['OS.time'], errors='coerce')
clin_clean['OS'] = pd.to_numeric(clin_clean['OS'], errors='coerce')
clin_clean = clin_clean.dropna()

clin_clean['gene_expr'] = clin_clean['sample'].map(expr_kirp)

# 6. Split by median
median_val = clin_clean['gene_expr'].median()
clin_clean['group'] = (clin_clean['gene_expr'] > median_val).map({True: 'High', False: 'Low'})

print(f"✅ High group: {sum(clin_clean['group'] == 'High')} patients")
print(f"✅ Low group: {sum(clin_clean['group'] == 'Low')} patients")

# 7. Kaplan-Meier plot
print("\n[5/5] Generating Kaplan-Meier plot...")
fig, (ax, ax_table) = plt.subplots(2, 1, figsize=(10, 8),
                                    gridspec_kw={'height_ratios': [3, 1]})

kmf_high = KaplanMeierFitter().fit(clin_clean[clin_clean['group'] == 'High']['OS.time'],
                                    clin_clean[clin_clean['group'] == 'High']['OS'],
                                    label='High ELFN1-AS1')
kmf_low = KaplanMeierFitter().fit(clin_clean[clin_clean['group'] == 'Low']['OS.time'],
                                   clin_clean[clin_clean['group'] == 'Low']['OS'],
                                   label='Low ELFN1-AS1')

kmf_high.plot_survival_function(ax=ax, color='red', ci_alpha=0.2)
kmf_low.plot_survival_function(ax=ax, color='blue', ci_alpha=0.2)

p_val = logrank_test(clin_clean[clin_clean['group'] == 'High']['OS.time'],
                      clin_clean[clin_clean['group'] == 'Low']['OS.time'],
                      event_observed_A=clin_clean[clin_clean['group'] == 'High']['OS'],
                      event_observed_B=clin_clean[clin_clean['group'] == 'Low']['OS']).p_value

ax.text(0.6, 0.1, f'p = {p_val:.4f}', transform=ax.transAxes, fontsize=14)
ax.set_title(f'Survival in {CANCER} by ELFN1-AS1 Expression')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Survival Probability')
ax.legend()
add_at_risk_counts(kmf_high, kmf_low, ax=ax_table, xticks=ax.get_xticks())
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/{CANCER}_survival.png", dpi=300)
plt.show()

print("\n" + "="*70)
print("✅ ANALYSIS COMPLETE")
print("="*70)
# =============================================================================
# Simplified and Corrected Kaplan-Meier Plot
# Author: Lubanah Younes
# =============================================================================

import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

# 1. Prepare data
high_time = clin_clean[clin_clean['group'] == 'High']['OS.time']
high_event = clin_clean[clin_clean['group'] == 'High']['OS']
low_time = clin_clean[clin_clean['group'] == 'Low']['OS.time']
low_event = clin_clean[clin_clean['group'] == 'Low']['OS']

# 2. Fit curves
kmf_high = KaplanMeierFitter().fit(high_time, high_event)
kmf_low = KaplanMeierFitter().fit(low_time, low_event)

# 3. Plot
plt.figure(figsize=(8,6))
kmf_high.plot_survival_function(color='red', label=f'High ELFN1-AS1 (n={len(high_time)})')
kmf_low.plot_survival_function(color='blue', label=f'Low ELFN1-AS1 (n={len(low_time)})')

# 4. Labels and title
plt.title('Survival in KIRP by ELFN1-AS1 Expression', fontweight='bold')
plt.xlabel('Time (days)')
plt.ylabel('Survival Probability')
plt.legend()

# 5. Add p-value
plt.text(0.6, 0.15, 'Log-rank p = 0.0018', 
         transform=plt.gca().transAxes, fontsize=12,
         bbox=dict(facecolor='white', alpha=0.8))

# 6. Save and show
plt.savefig(f"{OUTPUT_DIR}/KIRP_survival_corrected.png", dpi=300)
plt.show()
print("✅ Corrected plot saved.")

# =============================================================================
# Simple Cox Regression for ELFN1-AS1 in KIRP
# Author: Lubanah Younes
# =============================================================================

import pandas as pd
from lifelines import CoxPHFitter

print("\n" + "="*50)
print("SIMPLE COX REGRESSION (Univariate only)")
print("="*50)

# Prepare data directly from clin_clean
cox_simple = clin_clean[['OS.time', 'OS', 'gene_expr']].copy()
cox_simple.columns = ['time', 'event', 'gene_expr']
cox_simple = cox_simple.dropna()

print(f"Samples for univariate Cox: {cox_simple.shape[0]}")
print(f"Unique survival times: {cox_simple['time'].nunique()}")

if cox_simple['time'].nunique() > 1:
    cph = CoxPHFitter()
    cph.fit(cox_simple, duration_col='time', event_col='event')
    print("\n" + "="*30)
    print("Univariate Cox Results")
    print("="*30)
    print(cph.summary[['coef', 'exp(coef)', 'se(coef)', 'p']])
    
    # Save results
    cph.summary.to_csv(f"{OUTPUT_DIR}/KIRP_Cox_univariate.csv")
    print(f"\n✅ Results saved to {OUTPUT_DIR}/KIRP_Cox_univariate.csv")
else:
    print("⚠️ Cannot run Cox - all survival times are identical")

print("\n✅ Univariate analysis complete.")