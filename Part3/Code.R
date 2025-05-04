################################################################################ Start Load Seurat objects from Dataset 1: GSE144469
folder_path1 <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 1/6_GSE144469_seurat_objs"

rds_files1 <- list.files(folder_path1, pattern = "\\.rds$", full.names = TRUE)

for (i in seq_along(rds_files1)) {
  obj <- readRDS(rds_files1[i])
  assign(paste0("srob", i), obj)
}
################################################################################ End Load Seurat objects from Dataset 1: GSE144469



################################################################################ Start Load Seurat objects from Dataset 2: GSE206299
folder_path2 <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/5_GSE206299_seurat_objs"

rds_files2 <- list.files(folder_path2, pattern = "\\.rds$", full.names = TRUE)

offset <- length(rds_files1)  
for (i in seq_along(rds_files2)) {
  obj <- readRDS(rds_files2[i])
  assign(paste0("srob", i + offset), obj)
}
################################################################################ End Load Seurat objects from Dataset 2: GSE206299