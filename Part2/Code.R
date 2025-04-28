################################################################################ Start library
library(Seurat)
library(SeuratObject)
################################################################################ End library



################################################################################ Start CD4
h5ad_file <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/GSE206299_ircolitis-tissue-cd4.h5ad"
adata <- read_h5ad("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/GSE206299_ircolitis-tissue-cd4.h5ad")
mtx = adata$X
meta = adata$obs
genes = adata$var_names