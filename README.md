# ELFN1-AS1: A Novel Prognostic Biomarker in KIRP

## Overview
This repository contains the complete analysis of ELFN1-AS1 in Kidney Renal Papillary Cell Carcinoma (KIRP) using TCGA data. ELFN1-AS1 is identified as a novel prognostic biomarker with significant association with overall survival, cell cycle regulation, and immune infiltration.

## Key Findings
- Survival Analysis: Significant association with overall survival (Kaplan-Meier: p = 0.0018)
- Cox Regression: Independent prognostic factor (HR = 1.165, 95% CI: 1.047–1.296, p = 0.0049)
- Cell Cycle: Strong enrichment in G1/S-specific transcription (Odds Ratio = 15.47, p = 0.066)
- Immune Infiltration (TIMER3): 
  - Negative correlation with M2 macrophages (rho = -0.3578, p = 3.30e-09)
  - Positive correlation with MDSCs (rho = +0.2796, p = 5.12e-06)
  - Positive correlation with NK cells (rho = +0.2474, p = 5.87e-05)

## 📂 Repository Contents

### 📁 results_KIRP_ENSG00000238117/ – All output files
- ELFN1-AS1_Top7_Immune.csv – Top 7 immune cells (TIMER3 results)
- ELFN1-AS1_KIRP_Pathway_Results.xlsx – Pathway analysis (Reactome, KEGG)
- ELFN1-AS1_KIRP_Pathway_Results_1.jpg – Pathway visualization 1
- ELFN1-AS1_KIRP_Pathway_Results_2.jpg – Pathway visualization 2
- KIRP_Cox_univariate.csv – Univariate Cox regression results
- KIRP_survival.png – Kaplan-Meier survival plot
- KIRP_survival_corrected.png – Corrected survival plot
- ELFN1-AS1_top200_positive_genes.txt – Top 200 positively correlated genes
- ELFN1-AS1_top200_negative_genes.txt – Top 200 negatively correlated genes
- ELFN1-AS1_correlations.csv – Full correlation matrix

### 📜 Python Scripts
- KIRP_ELFN1-AS1_deep_analysis.py – Main analysis script (survival, Cox, correlations)
- extract_gene.py – Script to extract top correlated genes
- find_gene.py – Helper script to locate ELFN1-AS1 in expression data

### 📄 Supporting Files
- README.md – This file
- LICENSE – MIT-CR License
- requirements.txt – Python dependencies
- Survival_SupplementalTable_S1_20171025_xena_sp – TCGA clinical data (reference)
- tcga_RSEM_gene_tpm.gz – TCGA expression data (reference only)

## Requirements
- Python 3.7+
- Install dependencies: pip install -r requirements.txt

## Usage
1. Clone the repository.
2. Ensure data files are in the same directory.
3. Run the main analysis:
   `bash
   python KIRP_ELFN1-AS1_deep_analysis.py

د
## Citation
If you use this code or findings in your research, please cite:

Younes, L. (2026). ELFN1-AS1-KIRP: First release of ELFN1-AS1 KIRP analysis (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.19005670


## Author

LUBANAH YOUNES 081227

📧 lubanahyounes@gmail.com

## License

This project is licensed under the MIT-CR License – see the [LICENSE](LICENSE) file for details. Commercial use requires explicit permission from the author.
