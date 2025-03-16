library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(org.Hs.eg.db)
library(ggplot2)
library(GSEABase)
library(dplyr)
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library(RColorBrewer)

GDCprojects = getGDCprojects()

TCGAbiolinks:::getProjectSummary("TCGA-LIHC")

GDCprojects$id


query_LIHC = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  access = "open")

lihc_res = getResults(query_LIHC) # make results as table

GDCdownload(query_LIHC)
tcga_exp_LIHC=GDCprepare(query_LIHC)
saveRDS(tcga_exp_LIHC, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/LIHC_TCGA_data.rds")


compute_ECM_score <- function(se_obj, ecm_gene_set, exhaustion_gene_set) {
  # Ensure input is a SummarizedExperiment
  if (!inherits(se_obj, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object.")
  }
  
  # Extract metadata
  metadata <- colData(se_obj)
  filtered <- metadata$tumor_descriptor == "Primary" & 
    metadata$tissue_type == "Tumor"
  
  metadata_filtered <- metadata[filtered, ]
  
  expr_matrix <- SummarizedExperiment::assay(se_obj,'fpkm_uq_unstrand')[, filtered]
  # Extract Ensembl IDs from rownames of the expression matrix
  ensembl_ids <- rownames(expr_matrix)
  ensembl_ids_clean <- sub("\\..*", "", ensembl_ids)  # Remove version numbers
  # Update rownames in the expression matrix
  rownames(expr_matrix) <- ensembl_ids_clean
  
  gene_map <- AnnotationDbi::select(org.Hs.eg.db, keys = ensembl_ids_clean, 
                                    keytype = "ENSEMBL", columns = c("SYMBOL"))
  
  gene_map <- gene_map[!duplicated(gene_map$ENSEMBL), ]
  gene_map
  rownames(expr_matrix) <- gene_map$SYMBOL[match(rownames(expr_matrix), gene_map$ENSEMBL)]
  ## ECM
  ecm_expression <- expr_matrix[rownames(expr_matrix) %in% ecm_gene_set, ]
  missing_genes <- setdiff(ecm_gene_set, rownames(ecm_expression))
  if (length(missing_genes) > 0) {
    message("Warning: Some ECM genes are missing from dataset: ", paste(missing_genes, collapse = ", "))
  }
  ecm_stiffness_score <- colMeans(ecm_expression, na.rm = TRUE)
  
  # Create stiffness score dataframe
  ecm_score_df <- data.frame(
    Sample_ID = colnames(ecm_expression),
    Stiffness_Score = ecm_stiffness_score,
    Tumor_Type = metadata_filtered$project_id[match(colnames(ecm_expression), rownames(metadata_filtered))]
  )
  
  ## Exhaustion
  exhaustion_expression <- expr_matrix[rownames(expr_matrix) %in% exhaustion_gene_set, ]
  exhaustion_score <- colMeans(exhaustion_expression, na.rm = TRUE)
  
  exhaustion_df <- data.frame(
    Sample_ID = colnames(exhaustion_expression),
    Exhaustion_Score = exhaustion_score,
    Tumor_Type = metadata_filtered$project_id[match(colnames(exhaustion_expression), rownames(metadata_filtered))]
  )
  
  # Return both scores as a list
  return(list(stiffness_scores = ecm_score_df, exhaustion_scores = exhaustion_df))
}

#query_BRCA.2=query_BRCA
#tmp=query_BRCA.2$results[[1]]
#tmp=tmp[which(!duplicated(tmp$cases)),]
#query_BRCA.2$results[[1]]=tmp

DLBC_data <- readRDS("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/DLBC_TCGA_data.rds")
BRCA_data <- readRDS("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/BRCA_TCGA_data.rds")
ESCA_data <- readRDS("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/ESCA_TCGA_data.rds")
AML_data <- readRDS("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/LAML_TCGA_data.rds")
PACA_data <- readRDS("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/PAAD_TCGA_data.rds")
KIRC_data <- readRDS("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/KIRC_TCGA_data.rds")
THCA_data <- readRDS("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/THCA_TCGA_data.rds")
UCEC_data <- readRDS("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/UCEC_TCGA_data.rds")
OV_data <- readRDS("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/UCEC_TCGA_data.rds")
CHOL_data <- readRDS("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/CHOL_TCGA_data.rds")
LUAD_data <- readRDS("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/LUAD_TCGA_data.rds")
LIHC_data <- readRDS("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/LIHC_TCGA_data.rds")
COAD_data <- readRDS("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/COAD_TCGA_data.rds")
BLCA_data <- readRDS("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/BLCA_TCGA_data.rds")

genesets <- getGmt("EXTRACELLULAR_MATRIX_PART.v2024.1.Hs.gmt")
ecm_gene_list <- genesets[["EXTRACELLULAR_MATRIX_PART"]]@geneIds
# exhaustion_gene_list <- c("PDCD1", "LAG3", "TIGIT", "CTLA4", "HAVCR2", "TOX")
exhaustion_gene_list <- c("CD3D", "CD8A")
# exhaustion_gene_list <- c("CD274", "PDCD1LG2")
# exhaustion_gene_list <- c("ABHD12", "ACTA1", "ADGRV1", "ANGPT2", "ANKRD1", "ANKRD23", "ANO3", "AQP1", "ASIC2", "ASIC3",
#                               "ATAT1", "ATOH7", "ATP1A2", "ATP8A2", "ATR", "BACE1", "BAD", "BAG3", "BAK1", "BCL10",
#                               "BDKRB1", "BGLAP", "BMP6", "BNIP3", "BTG2", "CALB1", "CAPN2", "CASP1", "CASP2", "CASP5",
#                               "CASP8", "CASP8AP2", "CAV3", "CD40", "CDH2", "CHEK1", "CHI3L1", "CHRNA10", "CHRNA9", "CITED2",
#                               "CLCN6", "CNN2", "CNTNAP2", "COL11A1", "COL1A1", "COL6A1", "CRADD", "CSRP3", "CTNNB1", "CXCL10",
#                               "CXCL12", "CXCR4", "DAG1", "DCANP1", "DDR2", "DMD", "DRD2", "EDN1", "EGFR", "ENG",
#                               "ETV1", "F11R", "FADD", "FAS", "FGF2", "FOS", "FOSB", "FOSL1", "FYN", "GADD45A",
#                               "GAP43", "GATA4", "GCLC", "GDF5", "GPI", "GSN", "HABP4", "HPN", "HTR2A", "HTT",
#                               "IGF1R", "IGFBP2", "IHH", "IL13", "IL1B", "IL33", "IRF1", "ITGA2", "ITGAM", "ITGB3",
#                               "JUP", "KCNA1", "KCNA5", "KCNC1", "KCNJ2", "KCNK2", "KCNK4", "KCNQ1", "KCNQ3", "KIAA0319",
#                               "KIT", "KRT5", "LARGE1", "LHFPL5", "LRP11", "LTBR", "MAG", "MAP1B", "MAP2K4", "MAP3K1",
#                               "MAP3K14", "MAP3K2", "MAPK14", "MAPK3", "MAPK8", "MBD2", "MDK", "MEIS2", "MKKS", "MMP14",
#                               "MMP2", "MPO", "MTPN", "MYD88", "NEUROG1", "NFKB1", "NFKBIA", "NPPA", "NRXN1", "NRXN2",
#                               "NTRK1", "P2RX3", "P2RX7", "P2RY1", "PDE2A", "PDZD7", "PHF24", "PIEZO1", "PIEZO2", "PIK3CA",
#                               "PJVK", "PKD1", "PKD1L1", "PKD1L2", "PKD1L3", "PKD2", "PKD2L1", "PKD2L2", "PKDREJ", "PLEC",
#                               "POSTN", "PPL", "PSPH", "PTCH1", "PTGER4", "PTGS2", "PTK2", "PTK2B", "PTN", "RAF1",
#                               "RELA", "RETN", "RPS6KB1", "RYR2", "SCEL", "SCN11A", "SCN1A", "SCN9A", "SCX", "SERPINE2",
#                               "SHANK3", "SLC1A3", "SLC26A5", "SLC2A1", "SLC38A2", "SLC8A1", "SLC9A1", "SLITRK6", "SOST", "SOX9",
#                               "SRC", "STAT1", "STRA6", "STRBP", "STRC", "SUN1", "TACR1", "TCAP", "TGFB1", "THBS1",
#                               "TIFAB", "TLR3", "TLR4", "TLR5", "TLR7", "TLR8", "TMC1", "TMC2", "TMEM120A", "TMEM150C",
#                               "TMEM87A", "TNC", "TNF", "TNFRSF10A", "TNFRSF10B", "TNFRSF11A", "TNFRSF1A", "TNFRSF8", "TNFSF14", "TRPA1",
#                               "TRPV4", "TTN", "TUBA1A", "TXNIP", "UCN", "USP53", "WHRN", "WNT11", "XPA", "XPC" ) # Tumor malignancy Markers

print(exhaustion_gene_list %in% ecm_gene_list)
# Remove genes from exhaustion_gene_list if they exist in ecm_gene_list
exhaustion_gene_list <- exhaustion_gene_list[!(exhaustion_gene_list %in% ecm_gene_list)]
# Print the filtered gene list
print(exhaustion_gene_list)
print(ecm_gene_list)

DLBC_ecm_score <- compute_ECM_score(DLBC_data, ecm_gene_list,exhaustion_gene_list)
BRCA_ecm_score <- compute_ECM_score(BRCA_data, ecm_gene_list,exhaustion_gene_list)
ESCA_ecm_score <- compute_ECM_score(ESCA_data, ecm_gene_list,exhaustion_gene_list)
LAML_ecm_score <-  compute_ECM_score(AML_data, ecm_gene_list,exhaustion_gene_list)
PACA_ecm_score <- compute_ECM_score(PACA_data, ecm_gene_list,exhaustion_gene_list)
KIRC_ecm_score <- compute_ECM_score(KIRC_data, ecm_gene_list,exhaustion_gene_list)
THCA_ecm_score <- compute_ECM_score(THCA_data, ecm_gene_list,exhaustion_gene_list)
UCEC_ecm_score <- compute_ECM_score(UCEC_data, ecm_gene_list,exhaustion_gene_list)
OV_ecm_score <-   compute_ECM_score(OV_data, ecm_gene_list,exhaustion_gene_list)
CHOL_ecm_score <- compute_ECM_score(CHOL_data, ecm_gene_list,exhaustion_gene_list)
LUAD_ecm_score <- compute_ECM_score(LUAD_data, ecm_gene_list,exhaustion_gene_list)
LIHC_ecm_score <- compute_ECM_score(LIHC_data, ecm_gene_list,exhaustion_gene_list)
COAD_ecm_score <- compute_ECM_score(COAD_data, ecm_gene_list,exhaustion_gene_list)
BLCA_ecm_score <- compute_ECM_score(BLCA_data, ecm_gene_list,exhaustion_gene_list)

write.csv(DLBC_ecm_score$stiffness_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/DLBC_TCGA_stiffness_scores.csv", row.names = TRUE)
write.csv(DLBC_ecm_score$exhaustion_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/DLBC_TCGA_exhaustion_scores.csv", row.names = TRUE)
write.csv(BRCA_ecm_score$stiffness_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/BRCA_TCGA_stiffness_scores.csv", row.names = TRUE)
write.csv(BRCA_ecm_score$exhaustion_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/BRCA_TCGA_exhaustion_scores.csv", row.names = TRUE)
write.csv(ESCA_ecm_score$stiffness_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/ESCA_TCGA_stiffness_scores.csv", row.names = TRUE)
write.csv(ESCA_ecm_score$exhaustion_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/ESCA_TCGA_exhaustion_scores.csv", row.names = TRUE)
write.csv(LAML_ecm_score$stiffness_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/LAML_TCGA_stiffness_scores.csv", row.names = TRUE)
write.csv(LAML_ecm_score$exhaustion_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/LAML_TCGA_exhaustion_scores.csv", row.names = TRUE)
write.csv(PACA_ecm_score$stiffness_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/PACA_TCGA_stiffness_scores.csv", row.names = TRUE)
write.csv(PACA_ecm_score$exhaustion_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/PACA_TCGA_exhaustion_scores.csv", row.names = TRUE)
write.csv(KIRC_ecm_score$stiffness_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/KIRC_TCGA_stiffness_scores.csv", row.names = TRUE)
write.csv(KIRC_ecm_score$exhaustion_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/KIRC_TCGA_exhaustion_scores.csv", row.names = TRUE)
write.csv(THCA_ecm_score$stiffness_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/THCA_TCGA_stiffness_scores.csv", row.names = TRUE)
write.csv(THCA_ecm_score$exhaustion_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/THCA_TCGA_exhaustion_scores.csv", row.names = TRUE)
write.csv(UCEC_ecm_score$stiffness_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/UCEC_TCGA_stiffness_scores.csv", row.names = TRUE)
write.csv(UCEC_ecm_score$exhaustion_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/UCEC_TCGA_exhaustion_scores.csv", row.names = TRUE)
write.csv(OV_ecm_score$stiffness_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/OV_TCGA_stiffness_scores.csv", row.names = TRUE)
write.csv(OV_ecm_score$exhaustion_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/OV_TCGA_exhaustion_scores.csv", row.names = TRUE)
write.csv(CHOL_ecm_score$stiffness_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/CHOL_TCGA_stiffness_scores.csv", row.names = TRUE)
write.csv(CHOL_ecm_score$exhaustion_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/CHOL_TCGA_exhaustion_scores.csv", row.names = TRUE)
write.csv(LUAD_ecm_score$stiffness_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/LUAD_TCGA_stiffness_scores.csv", row.names = TRUE)
write.csv(LUAD_ecm_score$exhaustion_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/LUAD_TCGA_exhaustion_scores.csv", row.names = TRUE)
write.csv(LIHC_ecm_score$stiffness_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/LIHC_TCGA_stiffness_scores.csv", row.names = TRUE)
write.csv(LIHC_ecm_score$exhaustion_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/LIHC_TCGA_exhaustion_scores.csv", row.names = TRUE)
write.csv(COAD_ecm_score$stiffness_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/COAD_TCGA_stiffness_scores.csv", row.names = TRUE)
write.csv(COAD_ecm_score$exhaustion_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/COAD_TCGA_exhaustion_scores.csv", row.names = TRUE)
write.csv(BLCA_ecm_score$stiffness_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/BLCA_TCGA_stiffness_scores.csv", row.names = TRUE)
write.csv(BLCA_ecm_score$exhaustion_scores, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/BLCA_TCGA_exhaustion_scores.csv", row.names = TRUE)

stiffness_files <- list.files(path = "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/", 
                              pattern = "_TCGA_stiffness_scores.csv", full.names = TRUE)

exhaustion_files <- list.files(path = "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/", 
                               pattern = "_TCGA_exhaustion_scores.csv", full.names = TRUE)

# Extract the tumor type from the filename (e.g., "BLCA" from "BLCA_TCGA_stiffness_scores.csv")
stiffness_data <- bind_rows(lapply(stiffness_files, read.csv), .id = "File_Index")
exhaustion_data <- bind_rows(lapply(exhaustion_files, read.csv), .id = "File_Index")

# Extract Tumor Type from the filename
stiffness_data$Tumor_Type <- gsub("_TCGA_stiffness_scores.csv", "", basename(stiffness_files))[as.numeric(stiffness_data$File_Index)]
exhaustion_data$Tumor_Type <- gsub("_TCGA_exhaustion_scores.csv", "", basename(exhaustion_files))[as.numeric(exhaustion_data$File_Index)]

write.csv(stiffness_data, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/Stiffness_data_overall_data.csv", row.names = TRUE)
write.csv(exhaustion_data, "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/data_scores/Tcell_data_overall_data.csv", row.names = TRUE)
# combined_data <- merge(stiffness_data, exhaustion_data, by = c("Sample_ID", "Tumor_Type"))

# colnames(combined_data) <- c("Sample_ID", "Tumor_Type", "Stiffness_Score", "Exhaustion_Score")

# Compute mean stiffness score per tumor type
stiffness_rank <- stiffness_data %>%
  group_by(Tumor_Type) %>%
  summarise(Mean_Stiffness = mean(Stiffness_Score, na.rm = TRUE)) %>%
  arrange(Mean_Stiffness)

# Step 2: Convert Tumor_Type into a factor ranked by Mean_Stiffness
stiffness_data$Tumor_Type2 <- factor(stiffness_data$Tumor_Type, levels = stiffness_rank$Tumor_Type)

# stiffness_data$Tumor_Type <- factor(stiffness_data$Tumor_Type, levels = stiffness_rank$Tumor_Type)
# Set a professional color palette
# color_palette <- brewer.pal(n = min(15, length(unique(stiffness_data$Tumor_Type))), "Paired")

# Create Boxplot with Tumor Types Correctly Labeled on the X-Axis
fig1 <- ggplot(stiffness_data, aes(x = reorder(Tumor_Type, Stiffness_Score, fun = Mean), y = Stiffness_Score, fill = Tumor_Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  theme_minimal() +
  scale_fill_viridis_d(option = "plasma") +  # Automatically selects colors
  labs(title = "ECM Scores by Tumor Type",
       x = "Tumor Type (Ranked by Mean ECM Score)",
       y = "ECM Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))  # Adjust text size if needed

print(fig1)
ggsave("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/figures/Figure1_ECM_Stiffness_Boxplot.png", fig1,width = 10, height = 6, dpi = 600)


# Compute mean exhaustion score per tumor type
exhaustion_rank <- exhaustion_data %>%
  group_by(Tumor_Type) %>%
  summarise(Mean_Exhaustion = mean(Exhaustion_Score, na.rm = TRUE)) %>%
  arrange(Mean_Exhaustion)

# Plot Exhaustion scores by tumor type
fig2 <- ggplot(exhaustion_data, aes(x = reorder(Tumor_Type, Exhaustion_Score, fun = Mean), y = Exhaustion_Score, fill = Tumor_Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  theme_minimal() +
  scale_fill_viridis_d(option = "plasma") +  # Automatically selects colors
  labs(title = "Mechanosensation Scores by Tumor Type",
       x = "Tumor Type (Ranked by Mean Mechanosensation Score)",
       y = "Mechanosensation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))  # Adjust text size if needed

print(fig2)
ggsave("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/figures/Figure3_Tcell_Boxplot.png", fig2,width = 10, height = 6, dpi = 600)

# Compute mean scores per tumor type for correlation
correlation_data <- combined_data %>%
  group_by(Tumor_Type) %>%
  summarise(Mean_Stiffness = mean(Stiffness_Score, na.rm = TRUE),
            Mean_Exhaustion = mean(Exhaustion_Score, na.rm = TRUE))

# Correlation plot
fig3 <- ggplot(correlation_data, aes(x = Mean_Stiffness, y = Mean_Exhaustion, label = Tumor_Type)) +
  geom_point(size = 4, color = "blue") +
  geom_text(vjust = 1.5, hjust = 0.5, size = 4) +
  theme_minimal() +
  stat_smooth(method = "lm", color = "red", linetype = "dashed") + 
  labs(title = "Correlation Between Exhaustion Score and ECM Stiffness",
       x = "Mean ECM Stiffness Score",
       y = "Mean Exhaustion Score") +
  theme(axis.text = element_text(size = 12))

print(fig3)

ggsave("/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/figures/Figure1_ECM_Stiffness_Boxplot.png", fig1, width = 8, height = 6, dpi = 300)
ggsave("Figure2_Exhaustion_Score_Boxplot.png", fig2, width = 8, height = 6, dpi = 300)
ggsave("Figure3_Correlation_ECM_vs_Exhaustion.png", fig3, width = 7, height = 6, dpi = 300)

# rds_files <- list.files(path = "/oak/stanford/scg/lab_delitto/peter/ECM_Score_Bulk_TCGA/tcga_data/", pattern = "\\.rds$", full.names = TRUE)
# # Read all .rds files into a list
# tcga_list <- lapply(rds_files, readRDS)
# # write a for loop and then 
