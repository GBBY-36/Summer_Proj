# AML Transcriptome Preprocessing Pipeline

This repository contains a modular R preprocessing pipeline for Acute Myeloid Leukemia (AML) transcriptome data integration. It downloads, cleans, and merges RNA-Seq data from TCGA and GTEx (Discovery Cohort), and prepares the GSE13159 microarray dataset as an independent validation set.

## Project Structure

```text
├── .gitignore                   # Files and directories ignored by Git (e.g., data/ directories)
├── README.md                    # Project documentation (this file)
├── R_code.R                     # Synced master R script
├── R_code.Rmd                   # Synced master R Markdown notebook
├── group_meeting_report.docx    # Weekly report document (Word format)
└── scripts/                     # Modular preprocessing R scripts
    ├── 01_download_tcga_laml.R  # Downloader & cleaner for TCGA-LAML
    ├── 02_download_gtex.R       # Downloader & cleaner for GTEx Whole Blood
    ├── 03_download_gse13159.R   # Downloader & parser for GEO GSE13159
    └── 04_merge_discovery_cohort.R # Gene intersection & sample alignment tool
```

*Note: The raw data (`data_raw/`) and cleaned outputs (`data_clean/`) folders are excluded from Git repository to keep files lightweight and bypass GitHub file size limits.*

## Pipeline Steps

To reproduce the cleaned dataset from scratch, run the scripts in the terminal in the following sequence:

1. **Step 1: Download & Preprocess TCGA-LAML**
   ```bash
   Rscript scripts/01_download_tcga_laml.R
   ```
   *Downloads count matrices and clinical indices for 151 AML patients, cleans versioned Ensembl IDs, and scales TPM values to $\log_2(\text{TPM} + 1)$.*

2. **Step 2: Download & Preprocess GTEx Controls**
   ```bash
   Rscript scripts/02_download_gtex.R
   ```
   *Fetches GTEx V8 expression data, extracts 755 Whole Blood control samples, removes versioned Ensembl IDs, and scales to $\log_2(\text{TPM} + 1)$.*

3. **Step 3: Download & Preprocess Validation Set (GSE13159)**
   ```bash
   Rscript scripts/03_download_gse13159.R
   ```
   *Retrieves the MILE study dataset (2,096 samples), maps Affymetrix GPL570 probe IDs to Gene Symbols, aggregates duplicated genes (via data.table vectorization), and checks log scales.*

4. **Step 4: Merge Discovery Cohort**
   ```bash
   Rscript scripts/04_merge_discovery_cohort.R
   ```
   *Intersects shared genes (55,617 Ensembl genes) between TCGA and GTEx, column-binds the expression matrices, aligns and checks sample ID ordering dynamically.*

## Dependencies

All required package installations are automated at the beginning of each script. The dependencies are:
* **CRAN**: `data.table`, `dplyr`, `stringr`, `jsonlite`, `httr`, `plyr`
* **Bioconductor**: `TCGAbiolinks`, `SummarizedExperiment`, `GEOquery`, `Biobase`
