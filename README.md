# irAEs Project:
This repository contains a single-cell analysis pipeline for discovering biologically relevant insights in immune-related adverse events (irAE) disease.

---
---
---

## Overview:
The project is divided into four main parts:

1. **Part 1** – Analysis of the **GSE144469** dataset [NCBI GSE144469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144469) 
2. **Part 2** – Analysis of the **GSE206299** dataset [NCBI GSE206299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206299)
3. **Part 3** – Analysis of the **GSE206301** dataset [CPI-Colitis Single Cell and Spatial Atlas](https://simmonslab.shinyapps.io/CPI_COLITIS_DATA_PORTAL/)
4. **Part 4** – Integration of all three datasets

By analyzing each dataset individually and then integrating them, we aim to provide deeper insights into immune-related adverse event diseases.

---
---
---
## Data

The project's data summary is as follows:

### Dataset 1: GSE144469  
This dataset includes samples categorized into three cohorts:  
- **Normal controls**  
- **CPI-treated patients without colitis**  
- **CPI-treated patients with colitis**  

---

### Dataset 2: GSE206299  
This dataset includes samples from three distinct groups:  
- **Patients with irColitis**  
- **ICI-treated controls without irColitis**  
- **Healthy controls**  

---

### Dataset 3: GSE206301  
This dataset includes samples from five distinct groups:  
- **CPI Colitis**  
- **CPI Control**  
- **Healthy Control**  
- **UC Inflamed**  
- **UC Non-Inflamed**  

### Related Publications:
For more details on the data, please refer to the following publications:  
- [GSE144469](https://doi.org/10.1016/j.cell.2020.06.001)  
- [GSE206299](https://www.nature.com/articles/s41591-024-02895-x)  
- [GSE206301](https://www.nature.com/articles/s41591-024-02895-x)  

---
---
---

### Cell Type Focus:

- In the **first part** of the project, **CD3+ cells** were downloaded and analyzed from the GSE144469 dataset.  
- In the **second part**, **CD4+** and **CD8+ cells** were extracted and analyzed from the GSE206299 dataset.  
- In the **third part**, remaining cells from the GSE206301 dataset were selected based on upregulation of **CD3D**, **CD3E**, and **CD3G** genes, and downregulation of **LYZ**, **CD79A**, and **CD19** genes.
---
---
---

## Part 1 (GSE144469):
This part of the project consists of two sections. In the first section, as a warm-up, we draw a UMAP without data integration. In the second section, we integrate the data, draw the UMAP, and then identify gene markers and perform cell type annotation.

---

### First Step Analysis of Part 1:
As an initial analysis, we performed a **basic analysis** on the raw data without dataset integration. This step included:
- Standard Seurat preprocessing
- UMAP visualization of the merged object
- Identification of **marker genes** for each cluster

The related files, including plots and outputs, are saved in the **Merged_obj Data** folder.

---

### Second Step Analysis of Part 1:
In the next stage, we integrated the datasets using Seurat's integration workflow. The steps include:
- Performed quality control to ensure data integrity.
- Integrated the datasets to enable joint analysis.
- Generated UMAP plots for visualization of cellular distribution.
- Selected T cells
- Saved each sample separately as an individual Seurat object.
- Identified significant genes
- Annotated clusters based on identified marker genes, using both the provided Annotation Guide and supporting scientific literature to assign cell type identities.

---

### Annotation Guide
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

## Part 2: Adding the New Dataset (GSE206299):
In this section, we followed the same approach as in the second step of Part 1 of the project; however, this time, the analyses were performed on the GSE206299 dataset.

---
---
---

## Part 3: Integration of All 52 Samples Collected in Parts 1 and 2:
In this section, we applied the same analytical approach used in Parts 1 and 2 of the project. However, this time, the analyses were conducted on all 52 samples combined.