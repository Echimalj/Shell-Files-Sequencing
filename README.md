# human-aso-bulk-rnaseq

Bulk RNA-seq processing and differential expression analysis of human iPSC-derived neurons treated with tau antisense oligonucleotides (Tau ASO).

---

## Overview

This repository contains a reproducible pipeline for processing bulk RNA-seq data from human iPSC-derived glutamatergic neurons and performing downstream differential expression and overlap analysis.

The workflow includes:

- Alignment to the human genome (GRCh38) using STAR
- Gene-level quantification using RSEM
- Differential expression analysis using DESeq2
- Gene annotation (Ensembl → Symbol/Entrez)
- Cross-dataset overlap analysis with single-nucleus RNA-seq (snRNA-seq)

This pipeline was used to assess the concordance between tau-dependent transcriptional changes in human neuronal cultures and cell-type-specific signatures from mouse brain snRNA-seq.

---

## Data Source

Bulk RNA-seq data were obtained from:

- **GEO accession**: GSE204931

---

## Pipeline Summary

### 1. Genome Indexing (STAR)

```bash
sbatch Star_file.sh
```
### 2. Quantification (RSEM)

```bash
rsem-calculate-expression.sh
```
### 3. Count Matrix Generation

### 4. Differential Expression (DESeq2)

Performed in R
Performed in R:

Design: ~ condition
Contrast: Tau.ASO vs ASO.NT
Filtering: genes with total counts ≥ 10
Thresholds:
|log2FC| ≥ 0.58 (1.5-fold change)
adjusted p-value (BH) < 0.05

Outputs:

DESeq2_ASO_vs_Control.csv
DESeq2_results_annotated.csv

### 5. Gene Annotation
Gene IDs (Ensembl) are mapped to:

Gene symbols
Entrez IDs

### 6. Cross-Dataset Overlap Analysis
Differentially expressed genes (DEGs) from bulk RNA-seq are:

Mapped to mouse orthologs (via g:Profiler)
Compared to snRNA-seq-derived cell-type signatures

Overlap testing:

Fisher’s exact test (2 × 2 contingency tables)
Background universe:
all genes detected in snRNA-seq dataset
Multiple testing correction:
Benjamini-Hochberg (BH)


---
## Repository Structure

Shell-Files-Sequencing
├── scripts
│   ├── Star_file.sh                     # STAR genome indexing
│   ├── rsem-calculate-expression.sh     # RSEM quantification (Slurm array)
│   └── bulk_rnaseq_analysis.R           # DESeq2 + annotation workflow
│
├── README.md
└── .gitignore

---
### Notes
RSEM expected counts are rounded to integers before DESeq2 analysis.
Gene filtering removes lowly expressed genes to improve statistical power.
Annotation is based on Ensembl gene IDs.
Cross-species comparisons require ortholog mapping (not included in this repo).

