################################################################################ Start library
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(readxl)
library(future)
library(future.apply)
library(stringr)
library(cowplot)


if (!requireNamespace("MAST", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("MAST")
}
library(MAST)
################################################################################ End library



################################################################################ Start Load Seurat objects | Part 1
folder_path1 <- "C:/Esmaeil/Aim 2_irAEsProject/Backup/Part 1/3_The Seurat objects per sample"
rds_files1 <- list.files(folder_path1, pattern = "\\.rds$", full.names = TRUE)

for (i in seq_along(rds_files1)) {
  obj <- readRDS(rds_files1[i])
  assign(paste0("srob", i), obj)
}
################################################################################ End Load Seurat objects | Part 1


################################################################################ Start Load Seurat objects | Part 2
folder_path2 <- "C:/Esmaeil/Aim 2_irAEsProject/Backup/Part 2/5_The Seurat objects per sample"
rds_files2 <- list.files(folder_path2, pattern = "\\.rds$", full.names = TRUE)

offset2 <- length(rds_files1)
for (i in seq_along(rds_files2)) {
  obj <- readRDS(rds_files2[i])
  assign(paste0("srob", i + offset2), obj)
}
################################################################################ End Load Seurat objects | Part 2


################################################################################ Start Load Seurat objects | Part 3
folder_path3 <- "C:/Esmaeil/Aim 2_irAEsProject/Backup/Part 3/Without UC Non Inflamed/3_The Seurat objects per sample"
rds_files3 <- list.files(folder_path3, pattern = "\\.rds$", full.names = TRUE)

offset3 <- length(rds_files1) + length(rds_files2)

for (i in seq_along(rds_files3)) {
  obj <- readRDS(rds_files3[i])
  assign(paste0("srob", i + offset3), obj)
}
################################################################################ End Load Seurat objects | Part 3


################################################################################ Start Load Seurat objects | Part 4
folder_path4 <- "C:/Esmaeil/Aim 2_irAEsProject/Backup/Part 4/3_The Seurat objects per sample"
rds_files4 <- list.files(folder_path4, pattern = "\\.rds$", full.names = TRUE)

offset4 <- length(rds_files1) + length(rds_files2) + length(rds_files3)
for (i in seq_along(rds_files4)) {
  obj <- readRDS(rds_files4[i])
  assign(paste0("srob", i + offset4), obj)
}
################################################################################ End Load Seurat objects | Part 4



################################################################################ Start The fundamental steps before integration

merged_obj <- merge(srob1, y = list(srob2, srob3, srob4,srob5, srob6, srob7, srob8, srob9, srob10, srob11, srob12, srob13,
                                    srob14, srob15, srob16, srob17, srob18, srob19, srob20, srob21, srob22, srob23, srob24, srob25,
                                    srob26, srob27, srob28, srob29, srob30, srob31, srob32, srob33, srob34, srob35, srob36, srob37,
                                    srob38, srob39, srob40, srob41, srob42, srob43, srob44, srob45, srob46, srob47, srob48, srob49,
                                    srob50, srob51, srob52, srob53, srob54, srob55, srob56, srob57, srob58, srob59, srob60, srob61,
                                    srob62, srob63, srob64, srob65, srob66, srob67, srob68, srob69, srob70, srob71, srob72, srob73,
                                    srob74, srob75, srob76, srob77, srob78, srob79, srob80, srob81, srob82, srob83, srob84, srob85,
                                    srob86, srob87, srob88, srob89, srob90, srob91, srob92, srob93))
gc() 

rm(list = paste0("srob", 1:93))
gc() 

merged_obj=NormalizeData(merged_obj,normalization.method = "LogNormalize",scale.factor = 10000)
gc()

merged_obj=FindVariableFeatures(merged_obj,selection.method = "vst",nfeatures = 2000)
gc()

merged_obj = ScaleData(merged_obj)
gc()

merged_obj = RunPCA(merged_obj)
gc()

# 1
setwd("C:/Esmaeil/Aim 2_irAEsProject/Backup/Part 5/Without UC Non Inflamed/1_The Seurat object obtained after RunPCA and before IntegrateLayers")
saveRDS(file = "merged_obj",merged_obj)

################################################################################ End The fundamental steps before integration



################################################################################ Start integration with JointPCA

# ---------------------------------------------------- integration
library(future)
plan("sequential")
options(future.globals.maxSize = Inf)
options(future.fork.enable = FALSE)
options(future.rng.onMisuse = "ignore")

merged_obj <- IntegrateLayers(object = merged_obj,
                              method = JointPCAIntegration,
                              orig.reduction = "pca", 
                              new.reduction = "integrated.dr",
                              verbose = FALSE)

gc()


# 2
setwd("C:/Esmaeil/Aim 2_irAEsProject/Backup/Part 5/Without UC Non Inflamed/2_The Seurat object obtained after IntegrateLayers with JointPCAIntegration")
saveRDS(file = "merged_obj",merged_obj)


# ---------------------------------------------------- Post-integration Processing and Metadata Annotation
merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])


barcodes <- colnames(merged_obj)

Dataset <- str_split(barcodes, "/", simplify = TRUE)[,1]
Diagnosis_Patient <- str_split(barcodes, "/", simplify = TRUE)[,2]
Diagnosis <- str_split(Diagnosis_Patient, " ", simplify = TRUE)[,1]
Patient <- str_split(Diagnosis_Patient, " ", simplify = TRUE)[,2]

merged_obj@meta.data$Dataset <- Dataset
merged_obj@meta.data$Diagnosis <- Diagnosis
merged_obj@meta.data$Patient <- Patient


merged_obj1 = merged_obj


# ---------------------------------------------------- UMAP Visualization per Dataset for Integration QC
merged_obj1 = FindNeighbors(merged_obj1,dims = 1:30, reduction = "integrated.dr")
merged_obj1 = FindClusters(merged_obj1, resolution = 0.5)
merged_obj1 = RunUMAP(merged_obj1,dims = 1:30, reduction = "integrated.dr")


datasets <- unique(merged_obj1@meta.data$Dataset)

umap_list <- lapply(datasets, function(ds) {
  subset_obj <- subset(merged_obj1, subset = Dataset == ds)
  DimPlot(subset_obj, reduction = "umap") + ggtitle(paste("Dataset:", ds))
})

combined_plot <- plot_grid(plotlist = umap_list, ncol = 2)


ggsave("UMAP_by_Dataset_JointPCA.png", plot = combined_plot, width = 12, height = 10, dpi = 300)


# ---------------------------------------------------- Validate that cells belong to T-cell lineage
png(filename = "2_Feature Plot.png", width = 5000, height = 6000, units = "px", res = 600)
FeaturePlot(merged_obj1, features = c("CD3D", "CD3E", "CD3G", "LYZ", "CD79A", "CD19"))
dev.off()


# ---------------------------------------------------- Final UMAP visualization of integrated object
png(filename = "UMAP.png", width = 4000, height = 3000, units ="px", res = 600)
DimPlot(merged_obj1, label = TRUE)
dev.off()


################################################################################ End integration with JointPCA



################################################################################ Start integration with CCA

# ---------------------------------------------------- integration
plan("multicore", workers = 2)
options(future.globals.maxSize = 54 * 1024^3)

gc()
merged_obj <- IntegrateLayers(object = merged_obj,
                              method = CCAIntegration,
                              orig.reduction = "pca", 
                              new.reduction = "integrated.cca",
                              verbose = FALSE)

gc()

# 3
setwd("C:/Esmaeil/Aim 2_irAEsProject/Backup/Part 5/Without UC Non Inflamed/3_The Seurat object obtained after IntegrateLayers with CCAIntegration")
saveRDS(file = "merged_obj",merged_obj)


# ---------------------------------------------------- Post-integration Processing and Metadata Annotation
merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])


barcodes <- colnames(merged_obj)

Dataset <- str_split(barcodes, "/", simplify = TRUE)[,1]
Diagnosis_Patient <- str_split(barcodes, "/", simplify = TRUE)[,2]
Diagnosis <- str_split(Diagnosis_Patient, " ", simplify = TRUE)[,1]
Patient <- str_split(Diagnosis_Patient, " ", simplify = TRUE)[,2]

merged_obj@meta.data$Dataset <- Dataset
merged_obj@meta.data$Diagnosis <- Diagnosis
merged_obj@meta.data$Patient <- Patient


merged_obj1 = merged_obj


# ---------------------------------------------------- UMAP Visualization per Dataset for Integration QC
merged_obj = FindNeighbors(merged_obj,dims = 1:30, reduction = "integrated.cca")
merged_obj = FindClusters(merged_obj, resolution = 0.5)
merged_obj = RunUMAP(merged_obj,dims = 1:30, reduction = "integrated.cca")


datasets <- unique(merged_obj@meta.data$Dataset)

umap_list <- lapply(datasets, function(ds) {
  subset_obj <- subset(merged_obj, subset = Dataset == ds)
  DimPlot(subset_obj, reduction = "umap") + ggtitle(paste("Dataset:", ds))
})

combined_plot <- plot_grid(plotlist = umap_list, ncol = 2)

setwd("C:/Users/z5551702/OneDrive - UNSW/Aim 2_irAEsProject/irAEsProject/Part 5_Without UC Non Inflamed")
ggsave("UMAP_by_Dataset_CCA.png", plot = combined_plot, width = 12, height = 10, dpi = 300)


# ---------------------------------------------------- Validate that cells belong to T-cell lineage
png(filename = "2_Feature Plot_CCA.png", width = 5000, height = 6000, units = "px", res = 600)
FeaturePlot(merged_obj, features = c("CD3D", "CD3E", "CD3G", "LYZ", "CD79A", "CD19"))
dev.off()


# ---------------------------------------------------- Final UMAP visualization of integrated object
png(filename = "UMAP_CCA.png", width = 4000, height = 3000, units ="px", res = 600)
DimPlot(merged_obj, label = TRUE)
dev.off()


################################################################################ End integration with CCA

merged_obj2 = merged_obj1

################################################################################ Start Find Marker
merged_obj2 <- JoinLayers(merged_obj2)
setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part4")

Wilcoxmarkers = FindAllMarkers(merged_obj2,min.pct = 0.1 , logfc.threshold = 0.1, test.use = "wilcox")
write.csv(Wilcoxmarkers,file="WilcoxMarkers.csv")

MASTmarkers = FindAllMarkers(merged_obj2,min.pct = 0.1 , logfc.threshold = 0.1, test.use = "MAST")
write.csv(MASTmarkers,file="MASTMarkers.csv")
################################################################################ End Find Marker