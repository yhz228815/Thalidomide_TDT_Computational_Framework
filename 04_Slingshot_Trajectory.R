# =====================================================================
# Script 04: Pseudotime Trajectory Inference using Slingshot
# Description: Mapping the hematopoietic differentiation path
# =====================================================================

library(slingshot)
library(Seurat)
library(ggplot2)

# 1. Extract UMAP coordinates and cluster labels from Seurat object
umap_coords <- new_tdt_obj@reductions$umap@cell.embeddings
cell_labels <- new_tdt_obj@meta.data$lineages1

# 2. Fit trajectory, rooted at HSC/MPP
sds <- slingshot(umap_coords, clusterLabels = cell_labels, start.clus = "HSC/MPP")

# 3. Extract curve data for plotting
curves <- slingCurves(sds)
curve_data <- do.call(rbind, lapply(seq_along(curves), function(i) {
  data.frame(UMAP_1 = curves[[i]]$s[,1], UMAP_2 = curves[[i]]$s[,2], curve = as.character(i))
}))

# 4. Calculate Ribosomal Signature Score along Pseudotime
pseudotime <- slingPseudotime(sds)
pt_avg <- rowMeans(pseudotime, na.rm = TRUE)
expr_data <- GetAssayData(new_tdt_obj, assay = "RNA", layer = "data")
ribo_genes <- intersect(c("RPS3", "RPS9", "RPL34", "RPL39", "RPL11"), rownames(expr_data))
cell_scores <- colMeans(expr_data[ribo_genes, , drop = FALSE])

plot_df <- data.frame(Pseudotime = pt_avg, Ribosomal_Score = cell_scores)
plot_df <- plot_df[!is.na(plot_df$Pseudotime), ]