#!/usr/bin/env Rscript

# 03_dge_analysis.R
#
# This script performs differential gene expression analysis using DESeq2.
# It takes the raw count matrix from featureCounts and the sample metadata
# to identify genes that are differentially expressed between treatment conditions.
#
# It generates:
# 1. CSV files with DGE results for each comparison.
# 2. A PCA plot to visualize sample clustering.
# 3. A heatmap of the top differentially expressed genes.

# --- Load Libraries ---
# Note: You may need to install these packages first.
# install.packages(c("DESeq2", "tidyverse", "pheatmap", "ggrepel"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(ggrepel)
})

# --- Configuration ---
# Define project base directory relative to the script location
# This assumes you are running the script from the project root directory.
BASE_DIR <- getwd()

# Input files
COUNT_MATRIX_FILE <- file.path(BASE_DIR, "03_analysis/counts/raw_gene_counts.txt")
SAMPLE_SHEET_FILE <- file.path(BASE_DIR, "00_meta/sample_sheet.csv")

# Output directories
TABLES_DIR <- file.path(BASE_DIR, "03_analysis/tables")
FIGURES_DIR <- file.path(BASE_DIR, "03_analysis/figures")

# Create output directories if they don't exist
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Step 1: Load and Prepare Data ---
cat("Step 1: Loading and preparing data...\n")

# Load sample sheet
sample_sheet <- read_csv(SAMPLE_SHEET_FILE, show_col_types = FALSE) %>% 
  mutate(treatment = factor(treatment, levels = c("Untreated", "Dexamethasone", "Albuterol", "Albuterol_Dexamethasone")))

# Load count matrix from featureCounts
count_data <- read_tsv(COUNT_MATRIX_FILE, comment = "#", show_col_types = FALSE)

# --- Step 2: Clean and Format Count Matrix ---
# The column names from featureCounts are long file paths. We need to clean them
# to match the sample_id in the sample_sheet.
cat("Step 2: Cleaning and formatting the count matrix...\n")
clean_colnames <- colnames(count_data)[7:ncol(count_data)] %>%
  basename() %>%
  str_remove("_Aligned.sortedByCoord.out.bam") %>%
  str_remove("^SRR[0-9]+_") # <-- FIX: Remove the SRA run accession prefix

# Prepare the final count matrix with GeneID as row names
count_matrix <- count_data[, 7:ncol(count_data)] %>%
  as.matrix()
colnames(count_matrix) <- clean_colnames
rownames(count_matrix) <- count_data$Geneid

# --- Step 2.5: Robust Sanity Check ---
cat("Step 2.5: Verifying sample ID matching...\n")

# Check if all sample IDs from the sample sheet are present as columns in the count matrix
if (!all(sample_sheet$sample_id %in% colnames(count_matrix))) {
  missing_in_counts <- setdiff(sample_sheet$sample_id, colnames(count_matrix))
  stop(paste("Error: The following sample IDs from the sample sheet are NOT found in the count matrix columns:\n",
             paste(missing_in_counts, collapse="\n")))
}
# Check if all columns from the count matrix are present in the sample sheet
if (!all(colnames(count_matrix) %in% sample_sheet$sample_id)) {
  missing_in_sheet <- setdiff(colnames(count_matrix), sample_sheet$sample_id)
  stop(paste("Error: The following column names from the count matrix are NOT found in the sample sheet:\n",
             paste(missing_in_sheet, collapse="\n")))
}

# Ensure the column order in the count matrix matches the row order in the sample sheet
count_matrix <- count_matrix[, sample_sheet$sample_id]

# Final sanity check (this should now always pass)
stopifnot(all(colnames(count_matrix) == sample_sheet$sample_id))
cat("  - Sample IDs match perfectly.\n")


# --- Step 3: Create DESeq2 Object and Run Analysis ---
cat("Step 3: Creating DESeq2 object and running the analysis...\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_sheet,
  design = ~ treatment
)

# Pre-filter for low-count genes (e.g., keep genes with at least 10 reads total)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run the DESeq2 pipeline
dds <- DESeq(dds)

# --- Step 4: Generate and Save Diagnostic Plots ---
cat("Step 4: Generating diagnostic plots...\n")

# Variance Stabilizing Transformation for visualization
vsd <- vst(dds, blind = FALSE)

# Generate PCA plot
pca_plot <- plotPCA(vsd, intgroup = c("treatment", "cell_line")) +
  theme_bw() +
  ggtitle("PCA of Samples")
ggsave(file.path(FIGURES_DIR, "pca_plot.png"), pca_plot, width = 8, height = 6)
cat("  - PCA plot saved to", file.path(FIGURES_DIR, "pca_plot.png"), "\n")

# --- Step 5: Extract and Save DGE Results ---
cat("Step 5: Extracting and saving DGE results...\n")

# Define the comparisons of interest
comparisons <- list(
  Dex_vs_Untreated = c("treatment", "Dexamethasone", "Untreated"),
  Alb_vs_Untreated = c("treatment", "Albuterol", "Untreated"),
  AlbDex_vs_Untreated = c("treatment", "Albuterol_Dexamethasone", "Untreated")
)

# Loop through comparisons, get results, and save to CSV
for (comp_name in names(comparisons)) {
  res <- results(dds, contrast = comparisons[[comp_name]])
  res_df <- as.data.frame(res) %>%
    rownames_to_column(var = "gene_id") %>%
    arrange(padj)
  
  output_file <- file.path(TABLES_DIR, paste0(comp_name, "_dge_results.csv"))
  write_csv(res_df, output_file)
  cat("  - Results for", comp_name, "saved to", output_file, "\n")
}

# --- Step 6: Generate Heatmap of Top Genes ---
cat("Step 6: Generating a heatmap of top differentially expressed genes...\n")
# Using the Dex vs. Untreated comparison for this example
top_genes <- results(dds, contrast = comparisons$Dex_vs_Untreated) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  filter(!is.na(padj) & padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(30) # Select top 30 genes

# Get the normalized counts for these top genes
top_gene_counts <- assay(vsd)[top_genes$gene_id, ]

# Create the heatmap
pheatmap(top_gene_counts,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(dds)[, c("treatment", "cell_line")]),
         main = "Top 30 DE Genes (Dexamethasone vs. Untreated)",
         scale = "row",
         filename = file.path(FIGURES_DIR, "top30_DE_genes_heatmap.png"),
         width = 10, height = 12)
cat("  - Heatmap saved to", file.path(FIGURES_DIR, "top30_DE_genes_heatmap.png"), "\n")

cat("Differential expression analysis complete.\n")

