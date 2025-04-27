################################################################################ Start library
library(Seurat)
library(SeuratObject)


if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("zellkonverter")
library(zellkonverter)
################################################################################ End library



################################################################################ Start CD4
file_path <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/GSE206299_ircolitis-tissue-cd4.h5ad"
h5ad <- system.file("extdata", file_path ,package = "zellkonverter")
adata_CD4 = readH5AD(h5ad)