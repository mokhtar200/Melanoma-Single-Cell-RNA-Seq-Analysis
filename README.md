# Melanoma Single-Cell RNA-Seq Analysis

This repository contains a **full pipeline for analyzing melanoma single-cell RNA sequencing (scRNA-seq) data** using 10X Genomics data. The workflow includes **automatic cluster annotation with SingleR**, differential expression analysis, heatmaps, and functional enrichment (GO/KEGG).

---

## Project Overview
This project aims to **characterize the tumor microenvironment in melanoma at single-cell resolution**, focusing on immune and tumor cell populations. The main steps include:

1. Quality control and filtering of scRNA-seq data  
2. Dimensionality reduction (PCA, UMAP) and clustering  
3. Automatic cluster annotation using **SingleR**  
4. Differential expression analysis per cluster/cell type  
5. Functional enrichment (GO Biological Processes, KEGG pathways)
> **Note:** Keep raw data out of GitHub; use `.gitignore` for the `data/` folder.

---

## Requirements

- **R ≥ 4.2**
- R packages:
  - Seurat
  - SingleR
  - celldex
  - SingleCellExperiment
  - dplyr, ggplot2
  - clusterProfiler, enrichplot
  - org.Hs.eg.db

---

## Pipeline Overview

1. **Preprocessing & QC**
   - Filter low-quality cells
   - Normalize counts
   - Identify variable features

2. **Dimensionality Reduction & Clustering**
   - PCA → UMAP
   - Louvain clustering

3. **Automatic Cluster Annotation**
   - Use **SingleR** with Human Primary Cell Atlas reference
   - Assign biologically meaningful names to clusters

4. **Differential Expression Analysis**
   - Identify top marker genes per cluster/cell type
   - Visualize with heatmaps

5. **Functional Enrichment**
   - GO Biological Process
   - KEGG pathway analysis

---

