# Differential gene expression analysis using DESeq2

## 1. Build the DESeq2 dataset
```
R #if not already in R

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
dds <- DESeq(dds)
```
## 2. Get results, specify contrast explicitly for clarity - Ground Control vs. Flight
```
res <- results(dds, contrast = c("condition", "GC", "FLT"))
```
## 3. Convert results to a data frame for easier manipulation
```
res_df <- as.data.frame(res)
```
## 4. Add gene symbols
```
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Mm.eg.db")
library(AnnotationDbi)
library(org.Mm.eg.db)

res_df$symbol <- mapIds(org.Mm.eg.db,
                     keys = rownames(res),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")
```
## 5. Check for missing gene symbols - if it returns 0, that means no gene symbols are missing 
```
sum(is.na(res_df$symbol))
```
## 6. Order the results by adjusted significance/p-value
```
res_sig <- subset(res_df, padj < 0.05 & !is.na(padj))
```
## 7. Subset results for genes of interest: Ptgds and Cdkn2d
```
genes_of_interest <- c("Ptgds", "Cdkn2d")
subset_res <- subset(res_df, symbol %in% genes_of_interest)
print(subset_res)
```
## 8. Analyze the results

## p-values and adjusted p-values (padj) are given for each gene, as are the log2FoldChange values

## The adjusted p-value is used to correct for multiple testing - for example when you test thousands of genes - and reduces the false dicovery rate.  

## The log2FoldChange value measures the logarithmic ratio of gene expression levels between the two conditions - Flight vs. Ground Control. 

## In my analysis, positive log2FoldChange values = higher expression in Ground Control. Negative log2FoldChange values = higher expression in Flight.

# OPTIONAL:

## Look at the results for two additional genes to see if they were upregulated or downregulated as expected.

## GFAP is a classic marker for astrocyte activation (gliosis), which is a typical brain response to injury, stress, or inflammation - expected to be UPregulated in Flight

## BDNF supports neuron survival and synaptic plasticit - expected to be DOWNregulated in Flight.
```
genes_of_interest2 <- c("Gfap", "Bdnf")
subset_neuro_genes <- subset(res_df, symbol %in% genes_of_interest2)
print(subset_neuro_genes)
```
## Now it's time to visualize the results...

