# Germline Variant Calling Pipeline (GATK4 + Nextflow)

## Overview

This pipeline performs **germline short variant discovery (SNPs and INDELs)** from **paired-end human whole-genome sequencing (WGS) data** following the **GATK4 Best Practices** [workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)

It is a modular **Nextflow DSL2** reimplementation of the original bash [script](https://github.com/kpatel427/YouTubeTutorials/blob/main/variant_calling.sh), designed to be:
- Reproducible
- Scalable
- Containerised (Docker)
- Runnable locally or on AWS Batch

The workflow includes:
- Read QC
- Reference preparation
- Alignment
- Duplicate marking
- Base quality score recalibration (BQSR)
- Variant calling
- Variant selection (SNPs / INDELs)

---

## Pipeline Steps (Summary)

1. **FastQC** – Quality control of raw FASTQ files
2. **Reference preparation**
   - `samtools faidx`
   - `gatk CreateSequenceDictionary`
3. **Alignment** – BWA-MEM
4. **Mark duplicates** – GATK MarkDuplicatesSpark
5. **BQSR**
   - BaseRecalibrator
   - ApplyBQSR
6. **Metrics collection**
   - Alignment summary
   - Insert size metrics
7. **Variant calling** – GATK HaplotypeCaller
8. **Variant selection**
   - SNPs
   - INDELs

---

## Preliminary Steps

### 1. Install requirements

You need:

- **Nextflow** (≥ 22.x)
- **Docker** (for container execution)

Check:

```bash
nextflow -version
docker --version
```

### 2. Download example input data (from original script)

The original bash script uses paired-end reads from the 1000 Genomes Project.

```bash
wget -P data/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget -P data/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz
```

### 3. Download reference and known-sites files

```bash
# Reference genome
wget -P assets https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Known sites for BQSR
wget -P assets https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P assets https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
```

### 4. Use docker image to index ref file

```bash
# bwa index
docker run --rm \
  -v "$(pwd)/reference:/ref" \
  -w /ref \
  sktrinh12/bionfox-tools:1.0 \
  bwa index hg38.fa

# upload to s3 ref
aws s3 cp hg38.fa.fai s3://preludetx-strinh/nextflow/wgs-variant-calling-bioinformagician/reference/ --profile rchsandbox
aws s3 cp hg38.fa.amb s3://preludetx-strinh/nextflow/wgs-variant-calling-bioinformagician/reference/ --profile rchsandbox
aws s3 cp hg38.fa.ann s3://preludetx-strinh/nextflow/wgs-variant-calling-bioinformagician/reference/ --profile rchsandbox
aws s3 cp hg38.fa.bwt s3://preludetx-strinh/nextflow/wgs-variant-calling-bioinformagician/reference/ --profile rchsandbox
aws s3 cp hg38.fa.pac s3://preludetx-strinh/nextflow/wgs-variant-calling-bioinformagician/reference/ --profile rchsandbox
aws s3 cp hg38.fa.sa s3://preludetx-strinh/nextflow/wgs-variant-calling-bioinformagician/reference/ --profile rchsandbox
```

>Note: Reference indexing (`.fai`) and dictionary (`.dict`) are handled automatically by the pipeline.

Run BWA_INDEX() offline since it indexes the reference hg38.fa file and only need it once

```bash
docker run --rm \
  -v "$(pwd)/reference:/ref" \
  -w /ref \
  sktrinh12/bionfox-tools:1.0 \
  bwa index hg38.fa
[bwa_index] Pack FASTA... 20.94 sec
[bwa_index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=6418572210, availableWord=463634060
[BWTIncConstructFromPacked] 10 iterations done. 99999986 characters processed.
...
[BWTIncConstructFromPacked] 710 iterations done. 6418572210 characters processed.
[bwt_gen] Finished constructing BWT in 710 iterations.
[bwa_index] 5333.36 seconds elapse.
[bwa_index] Update BWT... 19.61 sec
[bwa_index] Pack forward-only FASTA... 14.12 sec
[bwa_index] Construct SA from BWT and Occ...
3950.69 sec
[main] Version: 0.7.19-r1273
[main] CMD: bwa index hg38.fa
[main] Real time: 9413.861 sec; CPU: 9338.763 sec
```

### Directory Layout

Recommended project structure:

```bash
.
├── assets
│  ├── hg38.fa
│  ├── Homo_sapiens_assembly38.dbsnp138.vcf
│  └── Homo_sapiens_assembly38.dbsnp138.vcf.idx
├── data
│  └── reads
│     ├── SRR062634_1.filt.fastq.gz
│     └── SRR062634_2.filt.fastq.gz
├── main.nf
├── modules
│  └── local
│     ├── apply_bqsr.nf
│     ├── base_recal.nf
│     ├── bwa_index.nf
│     ├── bwa_mem.nf
│     ├── collect_metrics.nf
│     ├── create_seq_dict.nf
│     ├── fastqc.nf
│     ├── haplotype_caller.nf
│     ├── mark_dups.nf
│     ├── samtool_faidx.nf
│     └── select_variants.nf
├── nextflow.config
└── results
   ├── aligned
   ├── metrics
   ├── qc
   └── variants
```

### Running the Pipeline

Local execution (Docker)

```bash
nextflow run main.nf -profile local
```

Resume if interrupted:

```bash
nextflow run main.nf -profile local -resume
```

### AWS Batch execution

Make sure:
* An AWS Batch queue exists

```bash

```
* An S3 bucket is configured as the Nextflow work directory

```bash
nextflow run main.nf -profile awsbatch
```

### Output

Results are written to the `results/` directory:

```bash
results/
├── qc/
├── aligned/
├── metrics/
└── variants/
    ├── raw_variants.vcf
    ├── raw_snps.vcf
    └── raw_indels.vcf
```
