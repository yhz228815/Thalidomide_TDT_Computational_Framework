# =====================================================================
# Script 05: Bulk RNA-seq Validation (GSE263306)
# Description: Comparing ribosomal signature score between TDT and Healthy
# =====================================================================

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

# 1. Load bulk data (Assuming pre-processed matrix is loaded)
# bulk_data <- read.table(gzfile("GSE263306_featureCounts_results_corrected_geo.txt.gz"), header = TRUE, row.names = 1)

# 2. Define core ribosomal genes
target_rps <- c("RPS3", "RPS9", "RPL34", "RPL39", "RPL11")
plot_data <- bulk_data[target_rps, , drop = FALSE]

# 3. Reshape and group (TDT vs Healthy)
plot_data$Gene <- rownames(plot_data)
plot_data_long <- pivot_longer(plot_data, 
                               cols = -c(Gene, Chr, Start, End, Strand, Length), 
                               names_to = "Sample", 
                               values_to = "Expression")
plot_data_long$Group <- ifelse(grepl("THAL", plot_data_long$Sample), "TDT", "Healthy")

# 4. Boxplot visualization
p_bulk <- ggplot(plot_data_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Healthy" = "#82B0D2", "TDT" = "#FA7F6F")) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE)