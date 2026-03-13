# ELFN1-AS1: A Novel Prognostic Biomarker in KIRP

This is the first comprehensive analysis of ELFN1-AS1 in KIRP, establishing it as a novel independent prognostic biomarker.

## Overview
This repository presents a comprehensive survival analysis of ELFN1-AS1 in Kidney Renal Papillary Cell Carcinoma (KIRP) using data from The Cancer Genome Atlas (TCGA). ELFN1-AS1 is a long non‑coding RNA (lncRNA) recently implicated in multiple cancers (colorectal, gastric, lung). This is the first study to identify ELFN1-AS1 as an independent prognostic biomarker in KIRP.

## Key Findings
- Significant association with overall survival (Kaplan–Meier: p = 0.0018)
- Independent prognostic factor in multivariate Cox regression  
  (Hazard Ratio = 1.165, 95% CI: 1.047–1.296, p = 0.0049)
- Analysis based on 321 KIRP samples with complete clinical and expression data from TCGA

## Repository Contents

### Core Scripts
- KIRP_ELFN1-AS1_deep_analysis.py – Main script for survival analysis (Kaplan–Meier, Cox regression, correlation analysis)
- KIRP_Pathway_Analysis.py – Pathway enrichment analysis for genes correlated with ELFN1-AS1
- find_gene.py – Helper script to locate ELFN1-AS1 in the expression file

### Results
- results_KIRP_ENSG00000238117/ – Folder containing all output files:
  - KIRP_Cox_univariate.csv – Univariate Cox regression results
  - KIRP_negative_pathways.csv – Pathway enrichment for negatively correlated genes
  - KIRP_positive_pathways.csv – Pathway enrichment for positively correlated genes
  - survival_ELFN1-AS1_KIRP.png – Kaplan–Meier survival plot

### Data Files
- Survival_SupplementalTable_S1_20171025_xena_sp – Clinical survival data from TCGA
- tcga_RSEM_gene_tpm.gz – Gene expression data (RSEM TPM) from TCGA

### Other
- requirements.txt – List of required Python packages
- LICENSE – MIT-CR license file
- README.md – This file

## Methodology
1. Data Loading  
   - Expression data: tcga_RSEM_gene_tpm.gz (gzipped)  
   - Clinical data: Survival_SupplementalTable_S1_20171025_xena_sp
2. Sample Filtering  
   - Retain only KIRP samples, remove those with missing survival data (OS.time or OS).
3. Gene Extraction  
   - Locate ELFN1-AS1 using its Ensembl ID ENSG00000238117 and extract expression values.
4. Patient Stratification  
   - Split patients into High and Low expression groups based on the median expression.
5. Survival Analysis  
   - Kaplan–Meier curves with log‑rank test (p = 0.0018).  
   - Cox proportional hazards regression (univariate) to obtain Hazard Ratio and significance (p = 0.0049).
6. Pathway Analysis (optional)  
   - Compute Pearson correlations between ELFN1-AS1 and all other genes.  
   - Submit top correlated genes to Enrichr for KEGG and GO enrichment analysis.

## Output Interpretation
- survival_ELFN1-AS1_KIRP.png – Survival plot showing the clear separation between High and Low groups.
- KIRP_Cox_univariate.csv – Contains:
  - coef : regression coefficient
  - exp(coef) : Hazard Ratio
  - p : p‑value (0.0049)
- KIRP_positive_pathways.csv / KIRP_negative_pathways.csv – Pathways significantly enriched among genes positively/negatively correlated with ELFN1-AS1.

## Requirements
- Python 3.7+
- Install dependencies:
 
  pip install -r requirements.txt

## Usage

1. Clone this repository:
  
   git clone https://github.com/Lubanah-Younes/ELFN1-AS1-KIRP.git
   cd ELFN1-AS1-KIRP
   
2. Ensure the required data files (tcga_RSEM_gene_tpm.gz and Survival_SupplementalTable_S1_20171025_xena_sp) are placed in the same directory as the scripts.
3. Run the main analysis:
  
   python KIRP_ELFN1-AS1_deep_analysis.py
   
4. (Optional) Run pathway analysis:
  
   python KIRP_Pathway_Analysis.py
   
## Citation

If you use this code or findings in your research, please cite:

· Younes, L. (2026). ELFN1-AS1: A Novel Prognostic Biomarker in KIRP. GitHub repository.
    https://github.com/Lubanah-Younes/ELFN1-AS1-KIRP

## Author

LUBANAH YOUNES 081227

📧 lubanahyounes@gmail.com

## License

This project is licensed under the MIT-CR License – see the [LICENSE](LICENSE) file for details. Commercial use requires explicit permission from the author.