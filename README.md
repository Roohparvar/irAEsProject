# irAEs Project:
This repository contains a single-cell analysis pipeline for discovering biologically relevant insights in immune-related adverse events (irAE) disease.

---
---
---

## Overview:
The project is divided into six main parts:

1. **Part 1** – Analysis of the **GSE144469** dataset [NCBI GSE144469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144469) 
2. **Part 2** – Analysis of the **GSE206299** dataset [NCBI GSE206299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206299)
3. **Part 3** – Analysis of the **GSE189040** dataset [CPI-Colitis Single Cell and Spatial Atlas](https://simmonslab.shinyapps.io/CPI_COLITIS_DATA_PORTAL/)
4. **Part 4** – Analysis of the **GSE253720** dataset [NCBI GSE253720](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253720)
5. **Part 5** – Analysis of the **GSE216329** dataset [NCBI GSE216329](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216329&utm_source=chatgpt.com)
6. **Part 6** – Integration of all five datasets

By analyzing each dataset individually and then integrating them, we aim to provide deeper insights into immune-related adverse event diseases.

---
---
---

## Data
The samples in this study were obtained through tissue biopsy collection.
The project's data summary is as follows:

### Dataset 1: GSE144469  
This dataset includes samples categorized into three cohorts:  
- **Normal controls**: Healthy individuals without immunotherapy treatment and no colitis (baseline control group).
- **CPI-treated patients without colitis**: Patients receiving Immune Checkpoint Inhibitors (CPI) treatment but without colitis.
- **CPI-treated patients with colitis**: Patients treated with CPI who developed colitis as an immune-related adverse event.  

---

### Dataset 2: GSE206299  
This dataset includes samples from three distinct groups:  
- **Patients with irColitis**: Healthy controls without any treatment and no colitis.
- **ICI-treated controls without irColitis**: Patients undergoing Immune Checkpoint Inhibitor (ICI) therapy without colitis. 
- **Healthy controls**: Patients on ICI therapy who developed immune-related colitis.

---

### Dataset 3: GSE189040  
This dataset includes samples from five distinct groups:  
- **CPI Colitis**: Patients treated with CPI who experienced irAE-associated colitis.
- **CPI Control**: Patients receiving CPI treatment but without colitis.
- **Healthy Control**: Healthy individuals with no treatment and no colitis. 
- **UC Inflamed**: Patients with active inflammation due to Ulcerative Colitis (UC). 
- **UC Non-Inflamed**: Patients with Ulcerative Colitis in remission or non-inflamed phase.

---

### Dataset 4: GSE253720  
This dataset includes samples from ... distinct groups:  
- ** ... ** : ... .
- ** ... ** : ... .
- ** ... ** : ... . 

---

### Dataset 5:   
This dataset includes samples from ... distinct groups:  
- ** ... ** : ... .
- ** ... ** : ... .
- ** ... ** : ... .

---

### Related Publications:
For more details on the data, please refer to the following publications:  
- [Part 1: GSE144469](https://doi.org/10.1016/j.cell.2020.06.001)  
- [Part 2: GSE206299](https://www.nature.com/articles/s41591-024-02895-x)  
- [Part 3: GSE189040](https://www.sciencedirect.com/science/article/pii/S153561082400134X?via%3Dihub#appsec2)  
- [Part 4: GSE253720](https://pubmed.ncbi.nlm.nih.gov/38642938/)  
- [Part 5: GSE216329](https://pubmed.ncbi.nlm.nih.gov/36513074/) 

---
---
---

## Cell Type Focus:

- In the **first part** of the project, **CD3+ cells** were downloaded and analyzed from the GSE144469 dataset.  
- In the **second part**, **CD4+** and **CD8+ cells** were extracted and analyzed from the GSE206299 dataset.  
- In the **third part**, cells were identified by high expression of **CD3D**, **CD3E**, and **CD3G** genes, and low expression of **LYZ**, **CD79A**, and **CD19** genes.

---
---
---

## Part 1: Analyzing Each Dataset Separately
For each dataset (article by article, dataset by dataset), we performed the following steps:

- Conducted quality control (QC) on individual samples  
- Integrated all samples within each dataset  
- Generated UMAP plots to visualize cellular distributions  
- Identified and selected T cell clusters  
- Saved each sample as a separate Seurat object for downstream analysis  

---
---
---

## Part 2: Integration of All Samples from the five articles
In this section, we applied the same analytical approach as in the previous parts. However, this time, the analysis was performed on all samples collectively by integrating the datasets obtained from the five referenced publications.

---
---
---

## Annotation Guide
Annotation is a critical step in the analysis process. The following guide is a starting point for identifying major immune cell types, though a thorough literature review is recommended for high-quality annotations.

#### Identifying Major Cell Types
- **CD4 cells**: Clusters expressing **"CD4"** and **"CD40LG"** markers.
- **CD8 cells**: Clusters expressing **"CD8A"** and **"CD8B"** markers.

#### CD4 Subclusters
##### Naïve CD4 T Cells
- **Upregulated genes**: `SELL`, `CCR7`
- **Downregulated genes**: `IL7R`

##### Regulatory T Cells (Treg)
- **Upregulated genes**: `CD4`, `CTLA4`, `FOXP3`, `TNFRSF4`, `TNFRSF18`, `TIGIT`, `IL2RA`
- **Downregulated genes**: Other markers not related to this subtype.

##### Follicular Helper T Cells (TFH)
- **Upregulated genes**: `CD4`, `ICOS`, `CD40LG`
- **Downregulated genes**: Other markers not related to this subtype.

##### Th1 Cells
- **Upregulated genes**: `CD4`, `CD40LG`, `IFNG`, `STAT1`, `STAT4`, `TBX21`
- **Downregulated genes**: Other markers not related to this subtype.

##### Th2 Cells
- **Upregulated genes**: `CD4`, `CD40LG`, `GATA3`, `STAT5`, `STAT6`
- **Downregulated genes**: Other markers not related to this subtype.

##### Th17 Cells
- **Upregulated genes**: `CD4`, `IL17A`, `RORA`, `CD40LG`
- **Downregulated genes**: Other markers not related to this subtype.

##### Effector Memory CD4 T Cells (CD4EM)
- **Upregulated genes**: `IL7R`, `LTB`, `CD40LG`, `CD4`
- **Downregulated genes**: `Sell`, `CCR7`

##### Central Memory CD4 T Cells (CD4CM)
- **Upregulated genes**: `IL7R`, `LTB`, `CD40LG`, `CD4`, `CCR7`
- **Downregulated genes**: `Sell`

##### Cytotoxic CD4 T Cells
- **Upregulated genes**: `CD4`, `CD4A`, `CD4B`, `GNLY`, `GZMA`, `GZMB`, `GZMH`, `GZMK`, `IFNG`, `PRF1`
- **Downregulated genes**: Other markers not related to this subtype.

#### CD8 Subclusters

##### Naïve CD8 T Cells
- **Upregulated genes**: `SELL`, `CCR7`
- **Downregulated genes**: Other markers not related to this subtype.

##### Cytotoxic CD8 T Cells
- **Upregulated genes**: `CD8A`, `CD8B`, `GNLY`, `GZMA`, `GZMB`, `GZMH`, `GZMK`, `IFNG`, `NKG7`, `PRF1`
- **Downregulated genes**: Other markers not related to this subtype.

#### Effector Memory CD8 T Cells (CD8EM)
- **Upregulated genes**: `IL7R`, `LTB`
- **Downregulated genes**: Specifically, `GNLY`, `GZMA`, `GZMB`, `GZMH`, and `GZMK`.

##### Central Memory CD8 T Cells (CD8CM)
- **Upregulated genes**: `IL7R`, `LTB`, `CD40LG`, `CD8A`, `CCR7`
- **Downregulated genes**: Other markers not related to this subtype.

---
---
---