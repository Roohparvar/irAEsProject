################################################################################ Start library
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(readxl)
################################################################################ End library



################################################################################ Start reading the RDS file and save each sample in a separate RDS file.
Sro_obj1 = LoadSeuratRds(file = "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part3/Data/CD3.single.cell.RDS")
setwd("C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part3/Data")

donor_info <- unique(Sro_obj1@meta.data[, c("Donor", "Type")])
type_matrix <- as.matrix(donor_info$Type)
colnames(type_matrix) <- "Type"

dir.create("Donor_Seurat_Objects", showWarnings = FALSE)
donor_info <- unique(Sro_obj1@meta.data[, c("Donor", "Type")])
group_counts <- list()

for (i in 1:nrow(donor_info)) {
  donor_id <- donor_info$Donor[i]
  donor_type <- gsub(" ", "_", donor_info$Type[i])
  
  if (!donor_type %in% names(group_counts)) {
    group_counts[[donor_type]] <- 1
  } else {
    group_counts[[donor_type]] <- group_counts[[donor_type]] + 1
  }
  
  suffix <- group_counts[[donor_type]]
  object_name <- paste0(donor_type, "_", suffix, "_New")
  
  donor_cells <- rownames(Sro_obj1@meta.data[Sro_obj1@meta.data$Donor == donor_id, ])
  seurat_obj <- subset(Sro_obj1, cells = donor_cells)
  
  assign(object_name, seurat_obj)
  
  saveRDS(seurat_obj, file = file.path("Donor_Seurat_Objects", paste0(object_name, ".rds")))
  
  print(paste("Saved:", object_name, ".rds"))
}
################################################################################ End reading the main RDS file and save each sample in a separate RDS file.



################################################################################ Start loading all RDS files and perform basic quality control on each Seurat object.
folder_path <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part3/Data/Donor_Seurat_Objects"
files <- list.files(path = folder_path, pattern = "\\.rds$", full.names = TRUE)
                

for (file in files) {
  base_name <- tools::file_path_sans_ext(basename(file))
  project_name <- base_name
  
  seurat_obj <- readRDS(file)
  
  if (!"nCount_RNA" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj[["nCount_RNA"]] <- Matrix::colSums(seurat_obj@assays$RNA@counts)
  }
  if (!"nFeature_RNA" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj[["nFeature_RNA"]] <- Matrix::colSums(seurat_obj@assays$RNA@counts > 0)
  }
  
  valid_cells <- WhichCells(seurat_obj, expression = nCount_RNA > 0 & nFeature_RNA > 0)
  seurat_obj <- subset(seurat_obj, cells = valid_cells)
  
  seurat_obj[["MTpercent"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  seurat_obj@project.name <- project_name
  
  assign(project_name, seurat_obj)
  
  
  output_base_path <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control"
  new_folder_path <- file.path(output_base_path, project_name)
  if (!dir.exists(new_folder_path)) {
    dir.create(new_folder_path)
  }
  
  cat("✅ Processed:", project_name, "with", length(valid_cells), "valid cells\n")
}
################################################################################ End loading all RDS files and perform basic quality control on each Seurat object.



################################################################################ Start CPI_Colitis_1_New

setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_1_New")
png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_1_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 2.5, color = 'red') + 
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 6000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_1_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 2.5, color = 'red') + 
  geom_vline(xintercept = 400, color = 'red') +
  geom_vline(xintercept = 2100, color = 'red')
dev.off()


CPI_Colitis_1_New=subset(CPI_Colitis_1_New,subset=nFeature_RNA>400 & nFeature_RNA<2100 & MTpercent<2.5 & nCount_RNA>700 & nCount_RNA<6000)

################################################################################ End CPI_Colitis_1_New



################################################################################ Start CPI_Colitis_2_New

setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_2_New")
png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_2_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 2.5, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 5000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_2_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 2.5, color = 'red') + 
  geom_vline(xintercept = 360, color = 'red') +
  geom_vline(xintercept = 1800, color = 'red')
dev.off()


CPI_Colitis_2_New=subset(CPI_Colitis_2_New,subset=nFeature_RNA>360 & nFeature_RNA<1800 & MTpercent<2.5 & nCount_RNA>600 & nCount_RNA<5000)

################################################################################ End CPI_Colitis_2_New



################################################################################ Start CPI_Colitis_3_New

setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_3_New")
       
png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_3_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 3, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 10000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_3_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 3, color = 'red') + 
  geom_vline(xintercept = 360, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


CPI_Colitis_3_New=subset(CPI_Colitis_3_New,subset=nFeature_RNA>360 & nFeature_RNA<3000 & MTpercent<3 & nCount_RNA>600 & nCount_RNA<10000)

################################################################################ End CPI_Colitis_3_New



################################################################################ Start CPI_Colitis_4_New

setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_4_New")
png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_4_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 2.5, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 3100, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_4_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 2.5, color = 'red') + 
  geom_vline(xintercept = 360, color = 'red') +
  geom_vline(xintercept = 1800, color = 'red')
dev.off()


CPI_Colitis_4_New=subset(CPI_Colitis_4_New,subset=nFeature_RNA>360 & nFeature_RNA<1800 & MTpercent<2.5 & nCount_RNA>600 & nCount_RNA<3100)

################################################################################ End CPI_Colitis_4_New



################################################################################ Start CPI_Colitis_5_New

setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_5_New")
png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_5_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 2.5, color = 'red') + 
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 2800, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_5_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 2.5, color = 'red') + 
  geom_vline(xintercept = 360, color = 'red') +
  geom_vline(xintercept = 1100, color = 'red')
dev.off()


CPI_Colitis_5_New=subset(CPI_Colitis_5_New,subset=nFeature_RNA>360 & nFeature_RNA<1100 & MTpercent<2.5 & nCount_RNA>700 & nCount_RNA<2800)

################################################################################ End CPI_Colitis_5_New





