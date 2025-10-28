# irAEs Project:
This repository presents a unified single-cell and TCR analysis pipeline focused on immune-related adverse events (irAEs). By analyzing multiple datasets of patients with and without checkpoint inhibitor–associated colitis, alongside healthy and ulcerative colitis controls, the project aims to uncover immune cell signatures, disease mechanisms, and potential biomarkers linked to irAEs.

---
---
---

## Overview:
The project is divided into five main parts:

1. **Part 1** – Analysis of the **GSE144469** dataset ([NCBI](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144469)) – Publication: [GSE144469](https://doi.org/10.1016/j.cell.2020.06.001)
2. **Part 2** – Analysis of the **GSE206299** dataset ([NCBI](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206299)) – Publication: [GSE206299](https://www.nature.com/articles/s41591-024-02895-x)
3. **Part 3** – Analysis of the **GSE189040** dataset ([CPI-Colitis Single Cell and Spatial Atlas](https://simmonslab.shinyapps.io/CPI_COLITIS_DATA_PORTAL/)) - Publication: [GSE189040](https://www.sciencedirect.com/science/article/pii/S153561082400134X)
4. **Part 4** – Analysis of the **GSE253720** dataset ([NCBI](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253720)) – Publication: [GSE253720](https://pubmed.ncbi.nlm.nih.gov/38642938/)
5. **Part 5** – Integration of all datasets, clustering, and cell type annotation.
6. **Part 6** – Downstream comparative and differential analysis.
7. **Part 7** – TCR clonotype and repertoire analysis.


---
---
---

## Data
The samples in this study were obtained through tissue biopsy collection.
The project's data summary is as follows:

### Dataset 1: GSE144469  
This dataset includes samples categorized into three cohorts:  
- **Healthy individuals with no treatment and no colitis**
- **Patients receiving Immune Checkpoint Inhibitors therapy  without colitisis**
- **Patients receiving Immune Checkpoint Inhibitors therapy who developed colitis**  

---

### Dataset 2: GSE206299  
This dataset includes samples from three distinct groups:  
- **Patients receiving Immune Checkpoint Inhibitors therapy who developed colitis**
- **Patients receiving Immune Checkpoint Inhibitors therapy  without colitis**
- **Healthy individuals with no treatment and no colitis**

---

### Dataset 3: GSE189040  
This dataset includes samples from five distinct groups:  
- **Patients receiving Immune Checkpoint Inhibitors therapy who developed colitis**
- **Patients receiving Immune Checkpoint Inhibitors therapy  without colitis**
- **Healthy individuals with no treatment and no colitis**
- **Patients with active inflammation due to Ulcerative Colitis (UC)**
- **Patients with Ulcerative Colitis in remission or non-inflamed phase**

---

### Dataset 4: GSE253720  
This dataset includes samples from ... distinct groups:  
- **Patients receiving Immune Checkpoint Inhibitors therapy who developed colitis**
- **Healthy individuals with no treatment and no colitis**
- **Patients with active inflammation due to Ulcerative Colitis (UC)**

---

### Data Grouping
The samples collected from the four datasets in this project are divided into **five main groups** based on treatment status, disease condition, and inflammation phase:

| Description | Short Name |
|------------|------------|
| Patients receiving Immune Checkpoint Inhibitors therapy who developed colitis | `CPI_Colitis` |
| Patients receiving Immune Checkpoint Inhibitors therapy  without colitis. | `CPI_Control` |
| Healthy individuals with no treatment and no colitis | `Healthy` |
| Patients with active inflammation due to Ulcerative Colitis (UC) | `UC_Inflamed` |
| Patients with Ulcerative Colitis in remission or non-inflamed phase | `UC_NonInflamed` |

---

### Related Publications:
For more details on the data, please refer to the following publications:  
- [Part 1: GSE144469](https://doi.org/10.1016/j.cell.2020.06.001)  
- [Part 2: GSE206299](https://www.nature.com/articles/s41591-024-02895-x)  
- [Part 3: GSE189040](https://www.sciencedirect.com/science/article/pii/S153561082400134X?via%3Dihub#appsec2)  
- [Part 4: GSE253720](https://pubmed.ncbi.nlm.nih.gov/38642938/)  

---
---
---

### Cell Type Focus:

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
In this section, we applied the same analytical approach as in the previous parts. However, this time, the analysis was performed on all samples collectively by integrating the datasets obtained from the five referenced publications.

---
---
---

