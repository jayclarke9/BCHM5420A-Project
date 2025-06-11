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
  
## Step 2: Download the quant.sf file
Download the **F1.quant.sf** file that we'll be using for this demo. If you remember from my talk, the `quant.sf` files are the output files that Salmon produces. They quantify the gene expression for RNA-Seq data. In my study I used 8 quant.sf files, but for this demo we'll just use one!
```r
download.file(
  "https://raw.githubusercontent.com/jayclarke9/BCHM5420A-Project/refs/heads/main/F1.quant.sf",
  destfile = "F1.quant.sf"
)

```
If the file downloaded correctly, you will see something like this:

`downloaded 5.5 MB`

## Step 3: Read this quant.sf file in R
```r
quant <- read.table("F1.quant.sf", header = TRUE, sep = "\t", quote = "")
```
## Step 4: Confirm the table looks good
View the first few rows of the table
```r
head(quant)
```
It should look like this:
                  Name Length EffectiveLength TPM NumReads
1 ENSMUST00000196221.2      9              10   0        0
2 ENSMUST00000179664.2     11              12   0        0
3 ENSMUST00000177564.2     16              17   0        0
4 ENSMUST00000178537.2     12              13   0        0
5 ENSMUST00000178862.2     14              15   0        0
6 ENSMUST00000179520.2     11              12   0        0
