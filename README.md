# Single-Cell Analysis Pipeline

This repository contains a single-cell analysis pipeline using the dataset **GSE144469**.

## Data

The overall design of the dataset includes samples from 22 patients across 3 different cohorts:

- **Normal control**: 8 people
- **CPI no colitis**: 6 people
- **CPI colitis**: 8 people

However, for the purposes of this project, we selected 3 samples from each of the 3 cohorts, meaning the project was conducted on 9 total samples:

- **Normal control**: 3 people
- **CPI no colitis**: 3 people
- **CPI colitis**: 3 people

It is important to note that we focused on **CD3+ cells** in this study.\
The data is available in the `Data` folder.

For more details on the methodology and findings, please refer to the published article:\
[https://www.nature.com/articles/s41591-024-02895-x](https://doi.org/10.1016/j.cell.2020.06.001)

To download the dataset, visit:\
[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144469)

## Focused Gene Sets

Based on previous studies, we identified well-known gene sets associated with **CD4+** and **CD8+ T cells**.\
We will focus on the expression of these genes during the `FindAllMarkers` step of the analysis.

### CD4+ T Cell Marker Genes

```
CD3D, CD4, CD40LG, TRAC, CD27, CD28, CD69, SELL, CCR7, LTB, IL7R, TCF7,
GATA3, TNF, JUN, IL32, IFNG, FOXP3, IL2RA, CTLA4, TIGIT, LAG3, PDCD1,
GZMA, GZMK, KLRB1, MKI67, TYMS, STMN1
```

### CD8+ T Cell Marker Genes

```
CD3D, CD8A, CD8B, TRAC, CD27, CD28, CD69, SELL, CCR7, LTB, IL7R, TCF7,
TBX21, IL32, PRF1, GZMA, GZMB, GZMK, IFNG, IFIT2, TNF, JUN, PDCD1,
LAG3, CTLA4, TIGIT, TRDC, TRGC2, CCL4, CCL5, GNLY, NKG7, TRAV1-2,
KLRB1, KLRD1, FCGR3A, FCER1G, NCAM1, MKI67, TYMS, STMN1
```

These gene sets will help us interpret the functional diversity of T cell subtypes and identify potential clusters with distinct immune profiles.

## Initial Analysis and Upcoming Steps

Please note that this project is ongoing, and updates are made regularly with new stages being added.
