---
# Single-Cell Analysis Pipeline
This repository contains a single-cell analysis pipeline using the dataset **GSE144469**.

## Data

The overall design of the dataset includes samples from 22 patients across 3 different cohorts:
- **Normal control**: 8 individuals  
- **CPI no colitis**: 6 individuals  
- **CPI colitis**: 8 individuals  

This study focuses specifically on **CD3+ cells**.  
The raw data can be found in the `Data` folder.

For more details on the methodology and findings, refer to the published article:  
[https://www.nature.com/articles/s41591-024-02895-x](https://doi.org/10.1016/j.cell.2020.06.001)

Dataset source (GSE144469):  
[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144469)

---

## First Step Analysis
As an initial and warm-up step, we performed a **basic analysis** on the raw data without dataset integration.  

This step includes:
- Standard Seurat preprocessing  
- UMAP visualization of the merged object  
- Identification of **marker genes** for each cluster

⚠️ All thresholds were chosen to support exploratory analysis and **can be modified** by users depending on their goals.

The related script is located in:  
`Without integration.R`  
Plots and outputs are saved in:  
`Merged_obj Data/`

Feel free to explore or extend this part of the analysis.

---

## Second Step Analysis
- The next stage includes integration of datasets using Seurat's integration workflow.
- A new UMAP is generated to compare cluster behavior across the integrated datasets, enabling improved cross-sample analysis and interpretation.
- Following integration and UMAP visualization, significant genes were identified using a defined threshold for adjusted p-value, average log fold change, and expression percentage differences.
- Based on the discovered marker genes and the guidance provided in the section 🧬 Annotation Guide—as well as through careful review of relevant scientific articles—annotation was performed to assign a cell type identity to each cluster.

---
## 🧬 Annotation Guide

**Annotation** is a critical phase that requires careful literature review and reading of relevant scientific articles. The guide below provides a helpful starting point, but **note**: this guide alone is not sufficient. A thorough literature review is essential to achieve high-quality annotations.

### 🔍 Identifying Major Cell Types

- If a cluster shows expression of **"CD4"** and **"CD40LG"** markers, it represents the **CD4 cell type**.
- If a cluster shows expression of **"CD8A"** and **"CD8B"** markers, it represents the **CD8 cell type**.

Both CD4 and CD8 T cells have several **subtypes**. Below, we first present **CD4 subclusters**, followed by **CD8 subclusters**.

#

### 🧪 CD4 Subclusters

#### • Naïve CD4 T Cells
- **Upregulated genes:** `SELL`, `CCR7`  
- **Downregulated genes:** `IL7R` (specifically), and all other markers not related to this subtype.

#### • Regulatory T Cells (Treg)
- **Upregulated genes:** `CD4`, `CTLA4`, `FOXP3`, `TNFRSF4`, `TNFRSF18`, `TIGIT`, `IL2RA`  
- **Downregulated genes:** All other markers not related to this subtype.

#### • Follicular Helper T Cells (TFH)
- **Upregulated genes:** `CD4`, `ICOS`, `CD40LG`  
- **Downregulated genes:** All other markers not related to this subtype.

#### • Th1 Cells
- **Upregulated genes:** `CD4`, `CD40LG`, `IFNG`, `STAT1`, `STAT4`, `TBX21`  
- **Downregulated genes:** All other markers not related to this subtype.

#### • Th2 Cells
- **Upregulated genes:** `CD4`, `CD40LG`, `GATA3`, `STAT5`, `STAT6`  
- **Downregulated genes:** All other markers not related to this subtype.

#### • Th17 Cells
- **Upregulated genes:** `CD4`, `IL17A`, `RORA`, `CD40LG`  
- **Downregulated genes:** All other markers not related to this subtype.

#### • Effector Memory CD4 T Cells (CD4EM)
- **Upregulated genes:** `IL7R`, `LTB`, `CD40LG`, `CD4`  
- **Downregulated genes:** `Sell, `CCR7` (specifically), and All other markers not related to this subtype.

#### • Central Memory CD4 T Cells (CD4CM)
- **Upregulated genes:** `IL7R`, `LTB`, `CD40LG`, `CD4`, `CCR7`  
- **Downregulated genes:** `Sell` (specifically), and all other markers not related to this subtype.

#### • Cytotoxic CD4 T Cells
- **Upregulated genes:** `CD4`, `CD4A`, `CD4B`, `GNLY`, `GZMA`, `GZMB`, `GZMH`, `GZMK`, `IFNG`, `PRF1`  
- **Downregulated genes:** All other markers not related to this subtype.

#

### 🧪 CD8 Subclusters

#### • Naïve CD8 T Cells
- **Upregulated genes:** `SELL`, `CCR7`  
- **Downregulated genes:** All other markers not related to this subtype.

#### • Cytotoxic CD8 T Cells
- **Upregulated genes:** `CD8A`, `CD8B`, `GNLY`, `GZMA`, `GZMB`, `GZMH`, `GZMK`, `IFNG`, `NKG7`, `PRF1`  
- **Downregulated genes:** All other markers not related to this subtype.

#### • Effector Memory CD8 T Cells (CD8EM)
- **Upregulated genes:** `IL7R`, `LTB`  
- **Downregulated genes:** Specifically, `GNLY`, `GZMA`, `GZMB`, `GZMH`, and `GZMK` should be downregulated. All other markers not related to this subtype should also be downregulated.

#### • Central Memory CD8 T Cells (CD8CM)
- **Upregulated genes:** `IL7R`, `LTB`, `CD40LG`, `CD8A`, `CCR7`  
- **Downregulated genes:** All other markers not related to this subtype should be downregulated.
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

## Upcoming Steps
This is an ongoing project.  
New stages and updates will be continuously added as the analysis progresses.

---
