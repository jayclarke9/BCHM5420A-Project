# Combine all the quant.sf files into one gene expression matrix

## 1. Creat a samples_metadata.csv file based on the existing salmon_output_* directories.
```
cd /home/paperspace/BCHM5420_Project/Data/Data

echo "sample_id,condition,path" > samples_metadata.csv

for dir in salmon_output_*; do
  sample=$(echo "$dir" | sed 's/salmon_output_//')

  if [[ "$sample" == *FLT* ]]; then
    condition="FLT"
  elif [[ "$sample" == *GC* ]]; then
    condition="GC"
  else
    condition="UNKNOWN"
  fi

  echo "$sample,$condition,$PWD/$dir/quant.sf" >> samples_metadata.csv
done
```
## 2. Check the samples_metadata.csv file to make sure it's correct
```
cat samples_metadata.csv
```
## 3. Download and unzip the mouse gtf file, which will be needed in step 6
```
cd /home/paperspace/BCHM5420_Project/Data/Data

wget ftp://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz
gunzip Mus_musculus.GRCm39.109.gtf.gz
```
## 4. Set up the R environment in VSCode
```
R
install.packages("BiocManager")
BiocManager::install(c("tximport", "readr", "DESeq2", "tximportData"))

library(tximport)
library(readr)
library(DESeq2)
```
## 5. Read the sample metadata
```
samples <- read.csv("samples_metadata.csv")
head(samples)
```
## 6. Create a tx2gene mapping file to translate data from transcript to gene
```
BiocManager::install("GenomicFeatures")
library(GenomicFeatures)

txdb <- makeTxDbFromGFF("Mus_musculus.GRCm39.109.gtf", format = "gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")

write.csv(tx2gene, "tx2gene.csv", row.names = FALSE)
```
## 7. Import counts using tximport
```
files <- samples$path
names(files) <- samples$sample_id

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
```
## Now that the data has been summarized as gene-level counts, proceed to differential gene expression analysis using DESeq2

