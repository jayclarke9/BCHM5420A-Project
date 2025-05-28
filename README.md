🚀🧠 ---- SpaceBrain ---- 🚀🧠

A Repository for Investigating Spaceflight-Induced Gene Expression Changes in the Rodent Brain

This repository contains a bioinformatics pipeline and downstream analysis for investigating gene expression changes in the hippocampus of rodents exposed to spaceflight. Specifically, the project focuses on the expression of Ptgds and Cdkn2d, two genes previously identified as significantly upregulated in spaceflight conditions. The analysis uses publicly available RNA-seq data from the NASA GeneLab database (Study OSD-612, GLDS-588).

This project uses the nf-core/rnaseq Nextflow pipeline to:
Perform quality control and trimming of raw RNA-seq data
Align reads and quantify gene expression
Generate a count matrix for downstream analysis
Assess differential expression of Ptgds and Cdkn2d using DESeq2

