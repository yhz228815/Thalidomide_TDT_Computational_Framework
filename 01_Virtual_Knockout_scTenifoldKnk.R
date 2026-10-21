# =====================================================================
# Script 01: Single-cell Virtual Knockout Analysis using scTenifoldKnk
# Description: Simulating the pharmacological inhibition of thalidomide targets
# =====================================================================

library(scTenifoldKnk)
library(Seurat)

# 1. Load the filtered TDT expression matrix
count_matrix_filtered <- readRDS("TDT_filtered_matrix_for_knk.rds")

# 2. Define core thalidomide targets for virtual knockout
target_genes <- c("CRBN", "BRD4", "IL1B")

# 3. Perform virtual knockout and save results
for (gene in target_genes) {
  cat("Running virtual knockout for:", gene, "\n")
  
  set.seed(123)
  result <- scTenifoldKnk(
    countMatrix = count_matrix_filtered,
    gKO = gene,
    qc_minLibSize = 500,
    nc_nNet = 5,
    nc_nCells = 500
  )
  
  # Extract significantly perturbed genes (FDR < 0.05)
  sig_genes <- subset(result$diffRegulation, p.adj < 0.05)
  saveRDS(result, file = paste0(gene, "_knockout_result.rds"))
}