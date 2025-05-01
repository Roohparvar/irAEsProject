# Single-Cell Analysis Pipeline

This repository contains a single-cell analysis pipeline for working with two datasets, **GSE144469** and **GSE206299**.

## Overview

In the first part of this project, we analyzed the **GSE144469** dataset. In the second part, we added the **GSE206299** dataset for a comprehensive analysis, combining both datasets to better understand immune cell populations across different conditions.

### Dataset 1: GSE144469

This dataset includes samples from 22 patients across 3 different cohorts:

- **Normal control**: 8 individuals
- **CPI no colitis**: 6 individuals
- **CPI colitis**: 8 individuals

### Dataset 2: GSE206299

The second dataset, added in the second part of the project, includes samples from 29 patients across 3 cohorts:

- **13 patients with irColitis** (9 on anti-PD-1 or anti-CTLA-4 monotherapy, and 4 on dual ICI therapy)
- **8 ICI-treated controls without irColitis**
- **8 healthy controls**

In the first part of the project, we focused specifically on **CD3+ cells**. In the second part, we focused on **CD4+** and **CD8+ cells**.


For more details on the methodology and findings, please refer to the following publications:
- [Nature Article on GSE144469](https://doi.org/10.1016/j.cell.2020.06.001)
- [Nature Article on GSE206299](https://www.nature.com/articles/s41591-024-02895-x)

Dataset Sources:
- **Part 1 (GSE144469)**: [NCBI GSE144469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144469)
- **Part 2 (GSE206299)**: [NCBI GSE206299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206299)

---

## First Step Analysis of Part 1

As an initial analysis, we performed a **basic analysis** on the raw data without dataset integration. This step included:

- Standard Seurat preprocessing
- UMAP visualization of the merged object
- Identification of **marker genes** for each cluster

⚠️ All thresholds were chosen to support exploratory analysis and **can be modified** based on user needs.

The related files, including plots and outputs, are saved in the **Merged_obj Data** folder.

---

## Second Step Analysis of Part 1

In the next stage, we integrated the datasets using Seurat's integration workflow. The steps include:

- Generation of a new UMAP to compare cluster behavior across the integrated datasets, enabling improved cross-sample analysis and interpretation.
- Identification of significant genes using thresholds for adjusted p-value, average log fold change, and expression percentage differences.
- Based on the identified marker genes, **annotation** was performed using both the provided **Annotation Guide** and a review of relevant scientific articles to assign cell type identities to each cluster.

---

## 🧬 Annotation Guide

Annotation is a critical step in the analysis process. The following guide is a starting point for identifying major immune cell types, though a thorough literature review is recommended for high-quality annotations.

### 🔍 Identifying Major Cell Types

- **CD4 cells**: Clusters expressing **"CD4"** and **"CD40LG"** markers.
- **CD8 cells**: Clusters expressing **"CD8A"** and **"CD8B"** markers.

### CD4 Subclusters

#### Naïve CD4 T Cells
- **Upregulated genes**: `SELL`, `CCR7`
- **Downregulated genes**: `IL7R`

#### Regulatory T Cells (Treg)
- **Upregulated genes**: `CD4`, `CTLA4`, `FOXP3`, `TNFRSF4`, `TNFRSF18`, `TIGIT`, `IL2RA`
- **Downregulated genes**: Other markers not related to this subtype.

#### Follicular Helper T Cells (TFH)
- **Upregulated genes**: `CD4`, `ICOS`, `CD40LG`
- **Downregulated genes**: Other markers not related to this subtype.

#### Th1 Cells
- **Upregulated genes**: `CD4`, `CD40LG`, `IFNG`, `STAT1`, `STAT4`, `TBX21`
- **Downregulated genes**: Other markers not related to this subtype.

#### Th2 Cells
- **Upregulated genes**: `CD4`, `CD40LG`, `GATA3`, `STAT5`, `STAT6`
- **Downregulated genes**: Other markers not related to this subtype.

#### Th17 Cells
- **Upregulated genes**: `CD4`, `IL17A`, `RORA`, `CD40LG`
- **Downregulated genes**: Other markers not related to this subtype.

#### Effector Memory CD4 T Cells (CD4EM)
- **Upregulated genes**: `IL7R`, `LTB`, `CD40LG`, `CD4`
- **Downregulated genes**: `Sell`, `CCR7`

#### Central Memory CD4 T Cells (CD4CM)
- **Upregulated genes**: `IL7R`, `LTB`, `CD40LG`, `CD4`, `CCR7`
- **Downregulated genes**: `Sell`

#### Cytotoxic CD4 T Cells
- **Upregulated genes**: `CD4`, `CD4A`, `CD4B`, `GNLY`, `GZMA`, `GZMB`, `GZMH`, `GZMK`, `IFNG`, `PRF1`
- **Downregulated genes**: Other markers not related to this subtype.

### CD8 Subclusters

#### Naïve CD8 T Cells
- **Upregulated genes**: `SELL`, `CCR7`
- **Downregulated genes**: Other markers not related to this subtype.

#### Cytotoxic CD8 T Cells
- **Upregulated genes**: `CD8A`, `CD8B`, `GNLY`, `GZMA`, `GZMB`, `GZMH`, `GZMK`, `IFNG`, `NKG7`, `PRF1`
- **Downregulated genes**: Other markers not related to this subtype.

#### Effector Memory CD8 T Cells (CD8EM)
- **Upregulated genes**: `IL7R`, `LTB`
- **Downregulated genes**: Specifically, `GNLY`, `GZMA`, `GZMB`, `GZMH`, and `GZMK`.

#### Central Memory CD8 T Cells (CD8CM)
- **Upregulated genes**: `IL7R`, `LTB`, `CD40LG`, `CD8A`, `CCR7`
- **Downregulated genes**: Other markers not related to this subtype.

---

## 🔍 Significant Gene Discovery

After identifying all marker genes, we applied the following thresholds to discover **significant genes**:

```
filtered_data <- subset(data, 
                        p_val < 0.05 & 
                        avg_log2FC > 0.1 & 
                        pct.1 > 0.1 & 
                        (pct.1 - pct.2) > 0.1)
```
---
# Part 2: Adding the New Dataset

In Part 2 of this project, we extended the analysis to include the **GSE206299** dataset, focusing on **CD4+** and **CD8+** cells in patients with **irColitis**. This dataset adds additional data points to the initial **GSE144469** dataset, allowing us to perform a more comprehensive immune cell analysis across different conditions.

---

## Dataset Overview: GSE206299

The **GSE206299** dataset includes samples from 29 patients, divided into the following cohorts:

- **irColitis** (treated with either anti-PD-1 or anti-CTLA-4 monotherapy or dual ICI therapy)
- **ICI-treated controls without irColitis**
- **healthy controls**

---

### Quality Control

As in the previous part of the analysis, we started with quality control steps. After integrating all data, we performed UMAP plotting. The key step in this part was to select only T cells.

### T Cell Selection

To identify T cells, we used specific marker genes: 
- **CD3D**, **CD3E**, **CD3G** for T cell identification
- **LYZ**, **CD78A** were not included in the T cell group.