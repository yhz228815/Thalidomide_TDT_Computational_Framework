# =====================================================================
# Script 03: Gene Ontology (GO) Enrichment Analysis
# Description: Identifying downstream pathways of perturbed genes
# =====================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# 1. Gene ID Conversion
gene_df <- bitr(all_sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 2. GO Biological Process Enrichment
ego <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# 3. Visualization (Dotplot & Barplot)
p_dot <- dotplot(ego, showCategory = 10)
p_bar <- barplot(ego, showCategory = 10)