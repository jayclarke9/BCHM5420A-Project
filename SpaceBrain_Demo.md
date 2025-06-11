# 🚀🧠 SpaceBrain Reproducibility Workshop 🧠🚀

In this short demo, we'll reproduce a single step from my project: extracting normalized expression (TPM) values from a Salmon `quant.sf` file using R.

In this demo, we’ll:
- Load a sample `quant.sf` file
- View its structure
- Extract the TPM value for *Ptgds*, a gene of interest in my project   
  

## Step 1: Install/Open R
- If you're **sure** you have R installed, open RStudio or run this in your terminal to open R:
```r
R
```
- If you're **unsure** if you have R installed, run this in your terminal to check. If you see something like R version 4.x.x, you’re good to go:
```r
R --version
```
- If you **need to install R**, you have a few options:

  - Activate our BCHM5420 class conda environment, which includes R:
  ```r
  curl -O https://raw.githubusercontent.com/uleth-advanced-bioinformatics/BCHM5420A-summer-2025/refs/heads/main/resources/advanced_bioinfo_env.yml
  conda env create -f advanced_bioinfo_env.yml
  conda activate uleth_advanced_bioinfo
  R
  ```
  - Install R onto your machine: https://cran.r-project.org
  
## Step 2: Download

    
