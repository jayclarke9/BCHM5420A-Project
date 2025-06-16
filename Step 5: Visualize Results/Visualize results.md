# Visualize the results using ggplot2 in R

# Barplot

## 1. In R, install the required libraries
```
install.packages("tidyr")
library(tidyr)
library(ggplot2)
library(DESeq2)
```
## 2. Set normalized counts
```
norm_counts <- counts(dds, normalized = TRUE)

genes_of_interest <- c("Ptgds", "Cdkn2d")
gene_ids <- rownames(res_df)[res_df$symbol %in% genes_of_interest]
```
## 3. Creat data frame for plotting
```
subset_counts <- norm_counts[gene_ids, ]
gene_symbols <- res_df[gene_ids, "symbol"]
df <- data.frame(
  gene = rep(gene_symbols, each = ncol(norm_counts)),
  sample = rep(colnames(norm_counts), times = length(gene_ids)),
  condition = rep(samples$condition, times = length(gene_ids)),
  count = as.vector(subset_counts)
)
```
# 4. Rename factor levels
```
df_long$condition <- factor(df_long$condition, 
                            levels = c("GC", "FLT"), 
                            labels = c("Ground Control", "Flight"))
```
# 5. Create and save plot as a JPEG
```
jpeg("gene_expression_barplot_final01.jpeg", width = 8, height = 6, units = "in", res = 300)

ggplot(df_long, aes(x = condition, y = count, fill = condition)) +
  stat_summary(geom = "bar", fun = mean, position = position_dodge(), width = 0.7) +
  stat_summary(geom = "errorbar", fun.data = mean_se, 
               position = position_dodge(.7), width = 0.25, linewidth = 0.5) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = c("Ground Control" = "#377EB8", "Flight" = "#E41A1C")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = NULL, y = "Normalized Counts")

dev.off()
```

# OPTIONAL: Volcano Plot showing all genes and their differential expression 

Ensure res_df has the necessary columns

If not already present, add gene symbols to the results data frame
```
res_df$gene <- res_df$symbol  # or modify as appropriate if symbols aren't yet present
```
Create a column for significance threshold
```
res_df$threshold <- with(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
```

Volcano plot
```
jpeg("volcano_plot_unique.jpeg", width = 8, height = 6, units = "in", res = 300)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = threshold), alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal(base_size = 14) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value",
    color = "Significant"
  ) +
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

dev.off()
```
