# Visualize the results using ggplot2 in R

# Barplot

## 1. In R, install the required libraries

install.packages("tidyr")
library(tidyr)
library(ggplot2)
library(DESeq2)

## 2. Set normalized counts
norm_counts <- counts(dds, normalized = TRUE)

genes_of_interest <- c("Ptgds", "Cdkn2d")
gene_ids <- rownames(res_df)[res_df$symbol %in% genes_of_interest]

## 3. Creat data frame for plotting
subset_counts <- norm_counts[gene_ids, ]
gene_symbols <- res_df[gene_ids, "symbol"]
df <- data.frame(
  gene = rep(gene_symbols, each = ncol(norm_counts)),
  sample = rep(colnames(norm_counts), times = length(gene_ids)),
  condition = rep(samples$condition, times = length(gene_ids)),
  count = as.vector(subset_counts)
)

# 4. Rename factor levels
df_long$condition <- factor(df_long$condition, 
                            levels = c("GC", "FLT"), 
                            labels = c("Ground Control", "Flight"))

# 5. Create and save plot as a JPEG
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


# Violin Plot showing all genes and their differential expression 





