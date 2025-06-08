################################################################################ Start library
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(readxl)
library(future)
library(future.apply)
################################################################################ End library



################################################################################ Start Load Seurat objects | Part 1
folder_path1 <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 1/5_GSE144469_seurat_objs"
rds_files1 <- list.files(folder_path1, pattern = "\\.rds$", full.names = TRUE)

for (i in seq_along(rds_files1)) {
  obj <- readRDS(rds_files1[i])
  assign(paste0("srob", i), obj)
}
################################################################################ End Load Seurat objects | Part 1


################################################################################ Start Load Seurat objects | Part 2
folder_path2 <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/5_GSE206299_seurat_objs"
rds_files2 <- list.files(folder_path2, pattern = "\\.rds$", full.names = TRUE)

offset2 <- length(rds_files1)
for (i in seq_along(rds_files2)) {
  obj <- readRDS(rds_files2[i])
  assign(paste0("srob", i + offset2), obj)
}
################################################################################ End Load Seurat objects | Part 2


################################################################################ Start Load Seurat objects | Part 3
folder_path3 <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part3/4_New_seurat_objs"
rds_files3 <- list.files(folder_path3, pattern = "\\.rds$", full.names = TRUE)

offset3 <- length(rds_files1) + length(rds_files2)
for (i in seq_along(rds_files3)) {
  obj <- readRDS(rds_files3[i])
  assign(paste0("srob", i + offset3), obj)
}
################################################################################ End Load Seurat objects | Part 3




################################################################################ Start integration
plan("multicore", workers = 2)
options(future.globals.maxSize = 54 * 1024^3)


merged_obj <- merge(srob1, y = list(srob2, srob3, srob4, 
                                    srob5, srob6, srob7, 
                                    srob8, srob9, srob10, 
                                    srob11, srob12, srob13,
                                    srob14, srob15, srob16,
                                    srob17, srob18, srob19,
                                    srob20, srob21, srob22,
                                    srob23, srob24, srob25,
                                    srob26, srob27, srob28,
                                    srob29, srob30, srob31,
                                    srob32, srob32, srob33,
                                    srob34, srob35, srob36,
                                    srob37, srob38, srob39,
                                    srob40, srob41, srob42,
                                    srob43, srob44, srob45,
                                    srob46, srob47, srob48,
                                    srob49, srob50, srob51,
                                    srob52, srob53, srob54,
                                    srob55, srob56, srob57,
                                    srob58, srob59, srob60,
                                    srob61, srob62, srob63,
                                    srob64, srob65, srob66,
                                    srob67, srob68, srob69,
                                    srob70, srob71, srob72,
                                    srob73, srob74, srob75,
                                    srob76, srob77, srob78,
                                    srob79, srob80, srob81,
                                    srob82, srob83))


rm(list = paste0("srob", 1:83))
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
setwd("C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part4/1_The Seurat object obtained after RunPCA and before IntegrateLayers")
saveRDS(file = "merged_obj",merged_obj)
# The Seurat object obtained after RunPCA and before IntegrateLayers
gc()

merged_obj <- IntegrateLayers(object = merged_obj,
                              method = CCAIntegration,
                              orig.reduction = "pca", 
                              new.reduction = "integrated.cca")

gc()
# 2
#setwd("C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part4/2_The Seurat object obtained after IntegrateLayers")
#saveRDS(file = "merged_obj",merged_obj)
# The Seurat object obtained after IntegrateLayers

merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])
################################################################################ End integration
