# 🧬 RNA-seq Pathway Activity Analysis using ssGSEA

**Author:** Oliver Abinader


## Overview

This repository provides a workflow to compute **pathway activity scores** from bulk RNA-seq data using **single-sample Gene Set Enrichment Analysis (ssGSEA)**.

Unlike differential expression analysis (gene-level), this approach summarizes expression into **biologically meaningful pathways**, enabling comparison of pathway activity across individual samples.


## Method Summary

* Input: TPM gene expression matrix
* Transformation: log(TPM + 1)
* Gene sets: Hallmark pathways (MSigDB, *Mus musculus*)
* Method: ssGSEA via `GSVA`
* Output:

  * Pathway activity scores per sample
  * Scaled (Z-score) pathway matrix
  * Heatmap visualization
  * Excel export


## Repository Structure

```
rnaseq-ssgsea-pathway-analysis/
│
├── scripts/
│   └── hallmark_ssgsea_analysis.Rmd
│
├── data/
│   └── metadata.tsv
│
├── results/
│   ├── Hallmark_ssGSEA.png
│   └── Hallmark_ssGSEA.xlsx
│
└── README.md
```


## 📥 Input Files

### 1. Expression Matrix

* File: `non-logged-scaled-TPM.allsamples.tsv`
* Format:

  * Rows = genes (GeneSymbol)
  * Columns = samples
  * Values = TPM

### 2. Metadata File

* File: `metadata.tsv`
* Required columns:

  * `sample`
  * `condition1` (e.g., group)
  * `condition2` (e.g., model/type)


## Workflow

### 1. Data Preprocessing

* Load TPM matrix
* Log-transform:

  ```
  log(TPM + 1)
  ```

### 2. Gene Set Retrieval

* Hallmark gene sets obtained using:

  * `msigdbr`
* Species: *Mus musculus*

### 3. ssGSEA Scoring

* Compute enrichment scores per:

  * sample
  * pathway
* Implemented using:

  ```r
  GSVA::gsva()
  ```

### 4. Scaling

* Z-score normalization per pathway:

  * highlights relative activation across samples

### 5. Visualization

* Heatmap generated using `pheatmap`
* Includes sample annotations:

  * Group
  * Model
    

## 📊 Outputs

| File                     | Description                  |
| ------------------------ | ---------------------------- |
| `Hallmark_ssGSEA.png`    | Heatmap of pathway activity  |
| `Hallmark_ssGSEA.xlsx`   | Raw and scaled ssGSEA scores |


## Interpretation

* Positive values → higher pathway activity
* Negative values → lower pathway activity


## ⚠️ Notes

* This is **not classical GSEA**
* ssGSEA operates at the **single-sample level**
* Results depend on:

  * normalization
  * gene set selection
  * data quality


## Dependencies

* `GSVA`
* `msigdbr`
* `dplyr`
* `tibble`
* `pheatmap`
* `openxlsx`


## 🧾 Citation

If you use this workflow, please cite:

* Hänzelmann et al., 2013 (GSVA method)
* MSigDB Hallmark gene sets


## Summary

This workflow transforms gene expression data into pathway-level activity scores, enabling biologically interpretable comparisons across samples.
