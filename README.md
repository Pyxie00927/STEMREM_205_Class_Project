# STEMREM_205_Class_Project
Single-cell Meta-analysis of T-cells Functional Phenotypes in Response to Extracellular Matrix Composition in Tumor Microenvironments 

Single-cell RNA seq of NSCLC
File_name = 'Analyses_NSCLC_T_cell_pancancer.ipynb'
Function: Preprocessing of single-cell RNA-seq dataset and analysis of exhaustion scores of T-cells in blood, adjacent normal tissue and tumor.

Bulk RNA-seq analysis
File_name = 'ECM_score_vPX.R'
Function: Score ECM genes, T-cell genes, mechanotransduction genes for 14 cancer types (fig. 2 of report)

Pre-processing of Mooney viscoelasticity dataset
File_name = 'Mooney_viscoelasticity_dataset_v3.ipynb'
Function: Preprocessing of single-cell RNA-seq dataset analysis and differential expression analysis between 3D collagen crosslinked (low viscoelasticity) with 3D collagen non-crosslinked (high viscoelasticity)

T-cell gene set: 
File_name = 'TcellGeneSet.Rmd'
Function: Use the DEGs from the Mooney viscoelasticity dataset to run GSEA and yield overlapping genes with response to mechanical stimulus pathway.

Pan-cancer Tcell analysis:
File_name = pan_canc_tcell_file_vfinal.ipynb
Function: Preprocessing of single-cell RNA-seq dataset and analysis of exhaustion and mechanotransduction scores of 10 different cancer types


