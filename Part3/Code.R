library(Seurat)
library(SeuratObject)
library(ggplot2)


################################################################################ Start Load Seurat objects from Dataset 1: GSE144469
folder_path1 <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Part 1_GSE144469_Seurat_Objs"

rds_files1 <- list.files(folder_path1, pattern = "\\.rds$", full.names = TRUE)

for (i in seq_along(rds_files1)) {
  obj <- readRDS(rds_files1[i])
  assign(paste0("srob", i), obj)
}
################################################################################ End Load Seurat objects from Dataset 1: GSE144469



################################################################################ Start Load Seurat objects from Dataset 2: GSE206299
folder_path2 <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Part 2_GSE206299_Seurat_Objs"

rds_files2 <- list.files(folder_path2, pattern = "\\.rds$", full.names = TRUE)

offset <- length(rds_files1)  
for (i in seq_along(rds_files2)) {
  obj <- readRDS(rds_files2[i])
  assign(paste0("srob", i + offset), obj)
}
################################################################################ End Load Seurat objects from Dataset 2: GSE206299


################################################################################ Start integration
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3")

merged_obj <- merge(srob1, y = list(srob2, srob3, srob4, 
                                    srob5, srob6, srob7, 
                                    srob8, srob9,srob10, 
                                    srob11, srob12,srob13,
                                    srob14, srob15,srob16,
                                    srob17, srob18,srob19,
                                    srob20, srob21,srob22,
                                    srob23, srob24,srob25,
                                    srob26, srob27,srob28,
                                    srob29, srob30,srob31,
                                    srob32, srob32,srob33,
                                    srob34, srob35,srob36,
                                    srob37, srob38,srob39,
                                    srob40, srob41,srob42,
                                    srob43, srob44,srob45,
                                    srob46, srob47,srob48,
                                    srob49, srob50,srob51,srob52))

gc()
merged_obj=NormalizeData(merged_obj,normalization.method = "LogNormalize",scale.factor = 10000)

gc()

merged_obj=FindVariableFeatures(merged_obj,selection.method = "vst",nfeatures = 2000)
#merged_obj=subset(merged_obj, features = VariableFeatures(merged_obj))
gc()

merged_obj = ScaleData(merged_obj)
gc()

merged_obj = RunPCA(merged_obj)
gc()

# 1
# saveRDS(file = "merged_obj",merged_obj)
# The Seurat object obtained after RunPCA and before IntegrateLayers
gc()
merged_obj <- IntegrateLayers(object = merged_obj,
                              method = CCAIntegration,
                              orig.reduction = "pca", 
                              new.reduction = "integrated.cca",
                              verbose = FALSE)

# 2
# saveRDS(file = "merged_obj",merged_obj)
# The Seurat object obtained after IntegrateLayers

merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])
################################################################################ End integration

merged_obj1 = merged_obj

################################################################################ Start UMAP

merged_obj1=FindNeighbors(merged_obj1,dims = 1:30,reduction = "integrated.cca")
merged_obj1=FindClusters(merged_obj1,resolution = 0.2)

# 3
# saveRDS(file = "merged_obj1",merged_obj1)
# The Seurat object obtained after FindNeighbors and FindClusters


merged_obj1=RunUMAP(merged_obj1,dims = 1:30, reduction = "integrated.cca")
png(filename = "UMAP1.png", width = 8000, height = 4000, units = "px", res = 1200)
DimPlot(merged_obj1, label = TRUE)
dev.off()

