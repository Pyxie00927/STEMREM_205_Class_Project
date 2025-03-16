# STEMREM 205 Class Project

## Single-cell Meta-analysis of T-cells Functional Phenotypes in Response to Extracellular Matrix Composition in Tumor Microenvironments

This repository contains code and analyses for the STEMREM 205 class project, which investigates T-cell functional phenotypes in response to extracellular matrix (ECM) composition within the tumor microenvironment. The project integrates single-cell and bulk RNA sequencing datasets across multiple cancer types to assess T-cell exhaustion and mechanotransduction responses.

---

## Project Files and Descriptions

### Single-cell RNA-seq of NSCLC
- **File:** `Analyses_NSCLC_T_cell_pancancer.ipynb`
- **Function:** Preprocessing of single-cell RNA-seq datasets and analysis of exhaustion scores of T-cells in blood, adjacent normal tissue, and tumor samples.

### Bulk RNA-seq Analysis
- **File:** `ECM_score_vPX.R`
- **Function:** Scores ECM genes, T-cell genes, and mechanotransduction genes for 14 cancer types (Figure 2 of the report).

### Pre-processing of Mooney Viscoelasticity Dataset
- **File:** `Mooney_viscoelasticity_dataset_v3.ipynb`
- **Function:** Preprocessing of single-cell RNA-seq dataset and differential expression analysis between 3D collagen crosslinked (low viscoelasticity) and 3D collagen non-crosslinked (high viscoelasticity) conditions.

### T-cell Gene Set Analysis
- **File:** `TcellGeneSet.Rmd`
- **Function:** Uses differentially expressed genes (DEGs) from the Mooney viscoelasticity dataset to run Gene Set Enrichment Analysis (GSEA) and identify overlapping genes with the response to mechanical stimulus pathway.

### Pan-cancer T-cell Analysis
- **File:** `pan_canc_tcell_file_vfinal.ipynb`
- **Function:** Preprocessing of single-cell RNA-seq datasets and analysis of exhaustion and mechanotransduction scores across 10 different cancer types.

### Correlation Score calculation
- **File:** `Correlation_scoring.ipynb`
- **Function:** Correlating the exhaustion and mechanotransduction scores from the pan-cancer dataset with the ECM scores from the TCGA dataset. 

---


