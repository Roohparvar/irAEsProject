# Single-Cell Analysis Pipeline
This repository contains a single-cell analysis pipeline using the dataset **GSE144469**.

## Data

The overall design of the dataset includes samples from 22 patients across 3 different cohorts:
- **Normal control**: 8 individuals  
- **CPI no colitis**: 6 individuals  
- **CPI colitis**: 8 individuals  

For the purposes of this project, we selected 3 samples from each cohort, resulting in a total of 9 samples:
- **Normal control**: 3 individuals  
- **CPI no colitis**: 3 individuals  
- **CPI colitis**: 3 individuals  

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
The next stage includes **integration of datasets** using Seurat's integration workflow.  
A new UMAP is generated to compare cluster behavior across the integrated datasets, enabling improved cross-sample analysis and interpretation.

---

## Upcoming Steps
This is an **ongoing project**.  
New stages and updates will be continuously added as the analysis progresses.