# 🚀🧠 ---- SpaceBrain ---- 🚀🧠

## A Repository for Investigating Spaceflight-Induced Gene Expression Changes in the Rodent Brain

This repository contains a bioinformatics pipeline and downstream analysis for investigating gene expression changes in the brain of rodents exposed to spaceflight. Specifically, the project focuses on the expression of Ptgds and Cdkn2d, two genes previously identified as significantly upregulated in spaceflight conditions. The analysis uses publicly available RNA-seq data from the NASA GeneLab database (Study OSD-612, GLDS-588).

This project uses the nf-core/rnaseq Nextflow pipeline to: perform quality control and trimming of raw RNA-seq data, align reads and quantify gene expression, generate a count matrix for downstream analysis, and assess differential expression of Ptgds and Cdkn2d using DESeq2.

## Workflow:

```mermaid
flowchart TD
    A("Start: Download FASTQ data from NASA GeneLab") --> B("Run nf-core/rnaseq pipeline with Nextflow")
    B --> C("Generate gene count matrix and QC reports")
    C --> D("Import count matrix into R")
    D --> E("Perform differential gene expression analysis using DESeq2")
    E --> F("Visualize Ptgds and Cdkn2d expression")
    F --> G("End: Insights into spaceflight-induced gene expression")

    style A fill:#C8E6C9
    style B fill:#C8E6C9
    style C fill:#C8E6C9
    style D fill:#C8E6C9
    style E fill:#C8E6C9
    style F fill:#C8E6C9
    style G fill:#C8E6C9
  ```
