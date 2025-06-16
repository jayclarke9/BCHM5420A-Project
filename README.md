# 🚀🧠 ---- SpaceBrain ---- 🧠🚀

## A Repository for Investigating Spaceflight-Induced Gene Expression Changes in the Rodent Brain

This repository contains a bioinformatics pipeline and downstream analysis for investigating gene expression changes in the brain of rodents exposed to spaceflight. Specifically, the project focuses on the expression of Ptgds and Cdkn2d, two genes previously identified as significantly upregulated in spaceflight conditions. The analysis uses publicly available RNA-seq data from the NASA GeneLab database (Study OSD-612, GLDS-588).

This project uses the nf-core/rnaseq Nextflow pipeline to: perform quality control and trimming of raw RNA-seq data, align reads and quantify gene expression, generate a count matrix for downstream analysis, and assess differential expression of Ptgds and Cdkn2d using DESeq2. 

An alternative quantification method using Salmon is also provided for users with limited computational resources.

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
## Quick Start Guide
### Option 1: Run full rnaseq Nextflow Pipeline
Ensure you have nextflow, conda (or docker), and nf-core/rnaseq installed.
```
nextflow run nf-core/rnaseq \
  --input samplesheet.csv \
  --genome GRCm38 \
  --outdir results/ \
  -profile conda
```
### Option 2: Lightweight Salmon quantification
For systems with limited memory or time constraints:
```
salmon quant -i transcripts_index \
  -l A \
  -1 sample_1.fastq.gz \
  -2 sample_2.fastq.gz \
  -p 4 \
  -o quants/sample_01
```
## Parameters and Definitions

### nf-core/rnaseq pipepline

| Parameter | Default Value | Description | Requirement |
|------------|---------|---------|--------------|
| `--input`   | `samplesheet.csv`  | CSV with sample info and FASTQ paths | Required |
| `--genome`  | None    | Reference ID from nf-core igenomes config (e.g. GRCm38 for mouse) | Required |
| `--outdir` | `./results` | Output directory | Optional |

### Salmon

| Parameter | Default Value | Description | Requirement |
|------------|---------|---------|--------------|
| `-i`   | N/A  | Path to Salmon transcriptome index | Required |
| `-l`  | A    | Library type (usually 'A' for autodetect') | Optional |
| `-1`, `-2` | N/A | Paired-end FASTQ files | Required |
| `-o` | `quants/` | Output directory for quantification results | Optional |
| `p` | `1` | Number of threads to use during quantification | Optional |

## Installation Instructions

### nf-core/rnaseq pipepline

1. Install Nextflow:
```
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/
```

2. Install Conda (for environments):
https://docs.conda.io/en/latest/miniconda.html

3. Launch pipeline:
```
nextflow run nf-core/rnaseq \
  --input samplesheet.csv \
  --genome GRCm38 \
  --outdir results/ \
  -profile conda
```
### Salmon

1. Install Salmon:
```
conda install -c bioconda salmon
```

2. Download and index the reference transcriptome:

You can download a mouse transcriptome FASTA file (e.g., from Ensembl or GENCODE), then build an index:
```
salmon index -t transcripts.fa -i transcripts_index
```

3. Run Salmon quantification:
```
salmon quant -i transcripts_index \
  -l A \
  -1 sample_1.fastq.gz \
  -2 sample_2.fastq.gz \
  -p 4 \
  -o quants/sample_01
```

## Outputs

### nf-core/rnaseq pipepline
```
results/
├── multiqc/
│   └── multiqc_report.html
├── star/
│   └── sample_01.bam
├── salmon/
│   └── quant.sf
├── counts/
│   └── gene_counts.tsv
└── logs/
    └── execution trace
```
| File | Description |
|------------|---------|
| `mutliqc_report.html` | Summary report of all QC tools |
| `gene_counts.tsv` | Count matrix for downstream expression analysis |
| `quant.sf` | Transcript-level quantification (TPM, read counts) |

### Salmon
```
quants/
├── sample_01/
│   ├── quant.sf
│   ├── cmd_info.json
│   └── lib_format_counts.json
```
| File | Description |
|------------|---------|
| `quant.sf` | Main quantification file used with 'tximport' in R |
| `lib_format_counts.json` | Count matrix for downstream expression analysis |
| `quant.sf` | Transcript-level quantification (TPM, read counts) |
