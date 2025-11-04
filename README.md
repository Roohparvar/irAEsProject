# irAEs Project:
This repository presents a unified single-cell and TCR analysis pipeline focused on immune-related adverse events (irAEs). By analyzing multiple datasets of patients with and without checkpoint inhibitor–associated colitis, alongside healthy and ulcerative colitis controls, the project aims to uncover immune cell signatures, disease mechanisms, and potential biomarkers linked to irAEs.

---
---
---

## Overview:
This project is organized into multiple parts, each focusing on the analysis of a different dataset, followed by integrated downstream and TCR analyses:

1. Analysis of the **GSE144469** dataset ([Access dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144469)) – and the [Publication](https://doi.org/10.1016/j.cell.2020.06.001)
2. Analysis of the **GSE206299** dataset ([Access dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206299)) – and the [Publication](https://www.nature.com/articles/s41591-024-02895-x)
3. Analysis of the **GSE189040** dataset ([Access dataset](https://simmonslab.shinyapps.io/CPI_COLITIS_DATA_PORTAL/)) – and the [Publication](https://www.sciencedirect.com/science/article/pii/S153561082400134X)
4. Analysis of the **GSE253720** dataset ([Access dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253720)) – and the [Publication](https://pubmed.ncbi.nlm.nih.gov/38642938/)
5. Integration of all datasets, clustering, and cell type annotation.
6. Downstream comparative and differential analysis.
7. TCR clonotype and repertoire analysis.


---
---
---

## Data
The samples in this study were obtained through tissue biopsy collection.
The project's data summary is as follows:

### Dataset 1: GSE144469  
This dataset includes samples categorized into three cohorts:  
- **Healthy individuals with no treatment and no colitis (8 Samples)**
- **Patients receiving Immune Checkpoint Inhibitors therapy  without colitis (6 Samples)**
- **Patients receiving Immune Checkpoint Inhibitors therapy who developed colitis (8 Samples)**  

##

### Dataset 2: GSE206299  
This dataset includes samples from three distinct groups:  
- **Healthy individuals with no treatment and no colitis(11 Samples)**
- **Patients receiving Immune Checkpoint Inhibitors therapy  without colitis(5 Samples)**
- **Patients receiving Immune Checkpoint Inhibitors therapy who developed colitis(14 Samples)**

##

### Dataset 3: GSE189040  
This dataset includes samples from five distinct groups:  
- **Healthy individuals with no treatment and no colitis (5 Samples)**
- **Patients receiving Immune Checkpoint Inhibitors therapy  without colitis (6 Samples)**
- **Patients receiving Immune Checkpoint Inhibitors therapy who developed colitis (7 Samples)**
- **Patients with active inflammation due to Ulcerative Colitis (9 Samples)**
- **Patients with Ulcerative Colitis in remission or non-inflamed phase (4 Samples)**

##

### Dataset 4: GSE253720  
This dataset includes samples from three distinct groups: 
- **Healthy individuals with no treatment and no colitis (3 Samples)** 
- **Patients receiving Immune Checkpoint Inhibitors therapy who developed colitis (9 Samples)**
- **Patients with active inflammation due to Ulcerative Colitis (2 Samples)**

##

### Data Grouping
The samples collected from the four datasets in this project are divided into **five main groups** based on treatment status, disease condition, and inflammation phase:

| Description | Short Name | Number of Samples | Number of Cells |
|------------|------------|-----------------|-----------------|
| Healthy individuals with no treatment and no colitis | `Healthy` | 27 | 59,845 |
| Patients receiving Immune Checkpoint Inhibitors therapy without colitis | `CPI_Control` | 17 | 41,398 |
| Patients receiving Immune Checkpoint Inhibitors therapy who developed colitis | `CPI_Colitis` | 38 | 92,088 |
| Patients with active inflammation due to Ulcerative Colitis (UC) | `UC_Inflamed` | 11 | 26,054 |
| Patients with Ulcerative Colitis in remission or non-inflamed phase | `UC_NonInflamed` | 4 | 14,845 |

##

### Cell Type Focus

- In the **first part** of the project, **CD3+ cells** were downloaded and analyzed from the GSE144469 dataset.  
- In the **second part**, **CD4+** and **CD8+ cells** were extracted and analyzed from the GSE206299 dataset.  
- In the **third and fourth parts**, cells were identified by high expression of **CD3D**, **CD3E**, and **CD3G** genes, and low expression of **LYZ**, **CD79A**, and **CD19** genes.

---
---
---

## Part 1: Analyzing Each Dataset Separately
For each dataset (article by article, dataset by dataset), we performed the following steps:

- Conducted quality control (QC) on individual samples. 
- Integrated all samples within each dataset.  
- Generated UMAP plots to visualize cellular distributions. 
- Identified and selected T cell clusters.
- Saved each sample as a separate Seurat object.

---
---
---

## Part 2: Integration of All Samples
In this section, we applied the same analytical approach as in the previous parts. However, this time, the analysis was performed on all samples collectively by integrating the datasets obtained from the all referenced datasets.

---
---
---

