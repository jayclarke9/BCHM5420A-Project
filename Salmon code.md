# RNA-Seq Quantification Pipeline Using Salmon

## 1. Activate Environment and Install Salmon
```bash
conda activate uleth_advanced_bioinfo
conda install -c bioconda salmon
```

## 2. Download mouse transcriptome and unzip it
```
get ftp://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
gunzip Mus_musculus.GRCm39.cdna.all.fa.gz
```

## 3. Create `sample_list.txt`
```bash
cd /home/paperspace/BCHM5420_Project/Data/Data

ls *_R1_raw.fastq.gz | sed 's/_R1_raw.fastq.gz//' > sample_list.txt

# Optional: view contents to check everything appears as expected (1 line per sample):
cat sample_list.txt
```

## 4. Build the Salmon index
```
salmon index -t Mus_musculus.GRCm39.cdna.all.fa -i mouse_salmon_index

```

## 5. Run Salmon Quantification in Parallel
```bash
cat sample_list.txt | parallel -j 4 '
  salmon quant -i /home/paperspace/mouse_salmon_index \
    -l A \
    -1 /home/paperspace/BCHM5420_Project/Data/Data/{}_R1_raw.fastq.gz \
    -2 /home/paperspace/BCHM5420_Project/Data/Data/{}_R2_raw.fastq.gz \
    -p 4 \
    -o /home/paperspace/BCHM5420_Project/Data/Data/salmon_output_{}
'
```
## Output
Will now have one folder per sample:

```
salmon_output_GLDS-588_rna-seq_RR10_BRN_FLT_WT_F1/
 └── quant.sf
```

Each `quant.sf` contains transcript-level quantifications (TPM, estimated counts, etc.).

---

## Ready for Downstream Analysis
Next steps:
- Aggregating with **tximport** in R
- Summarizing to gene-level
- Running differential expression analysis using **DESeq2**
