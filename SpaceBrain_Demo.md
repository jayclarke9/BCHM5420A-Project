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
Download the **F1.quant.sf** file that we'll be using for this demo. If you remember from my talk, the `quant.sf` files are the output files that Salmon produces. They quantify the gene expression from RNA-Seq data. In my study I used 8 quant.sf files, but for this demo we'll keep it simple and just use one!
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

![table example](https://github.com/user-attachments/assets/286a1c21-a419-48e7-b137-a41cfccc4671)

## Step 5: Find the TPM value for the *Ptgds* gene
As you can see from our table, the `Name` column doesn't have recognizable gene names, it has the Ensembl transcript ID for each gene (e.g. 'ENSMUST00000196221.2').

At this point in my study, I used a `tx2gene` file to map these transcript IDs to gene names. But we probably don't have time for that! 

Let's use the Ensembl IDs instead. A quick look at [Ensembl online](https://www.ensembl.org/Mus_musculus/Transcript/Summary?db=core;g=ENSMUSG00000015090;r=2:25356721-25359854;t=ENSMUST00000015234) tells us that the Ensembl ID for the Ptgds gene in mice is `ENSMUST00000015234`.

Now let's find the **TPM value** for that gene ID. TPM stands for **Transcripts Per Million** and is the  relative abundance measure that we use for downstream analysis.

To find the TPM value for our *Ptgds* gene, use this code:
```r

