Attempted rnaseq pipeline command with skipped processes and a custom config file:

```
nextflow run nf-core/rnaseq -profile docker --input ~/BCHM5420_Project/Data/Data/samplesheet.csv --outdir ~/BCHM5420_Project/Results --genome GRCm38 -c custom.config --skip_fastqc --skip_multiqc --skip_trimming  --aligner star_salmon
```
