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



################################################################################ Start CPI_Colitis_6_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_6_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_6_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 1450, color = 'red') +
  geom_vline(xintercept = 6500, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_6_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 2400, color = 'red')
dev.off()


CPI_Colitis_6_New=subset(CPI_Colitis_6_New,subset=nFeature_RNA>600 & nFeature_RNA<2400 & MTpercent<5 & nCount_RNA>1450 & nCount_RNA<6500)
################################################################################ End CPI_Colitis_6_New



################################################################################ Start CPI_Colitis_7_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_7_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_7_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 1100, color = 'red') +
  geom_vline(xintercept = 4500, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_7_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 500, color = 'red') +
  geom_vline(xintercept = 2000, color = 'red')
dev.off()


CPI_Colitis_7_New=subset(CPI_Colitis_7_New,subset=nFeature_RNA>500 & nFeature_RNA<2000 & MTpercent<5 & nCount_RNA>1100 & nCount_RNA<4500)
################################################################################ End CPI_Colitis_7_New



################################################################################ Start CPI_Colitis_8_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_8_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_8_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 1100, color = 'red') +
  geom_vline(xintercept = 4500, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_8_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 500, color = 'red') +
  geom_vline(xintercept = 2000, color = 'red')
dev.off()


# CPI_Colitis_8_New=subset(CPI_Colitis_8_New,subset=nFeature_RNA>500 & nFeature_RNA<2000 & MTpercent<5 & nCount_RNA>1100 & nCount_RNA<4500)
################################################################################ End CPI_Colitis_8_New



################################################################################ Start CPI_Colitis_9_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_9_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_9_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 1100, color = 'red') +
  geom_vline(xintercept = 4500, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_9_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 500, color = 'red') +
  geom_vline(xintercept = 2000, color = 'red')
dev.off()


# CPI_Colitis_9_New=subset(CPI_Colitis_9_New,subset=nFeature_RNA>500 & nFeature_RNA<2000 & MTpercent<5 & nCount_RNA>1100 & nCount_RNA<4500)
################################################################################ End CPI_Colitis_9_New



################################################################################ Start CPI_Colitis_10_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_10_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_10_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 4500, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_10_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 900, color = 'red') +
  geom_vline(xintercept = 1700, color = 'red')
dev.off()


CPI_Colitis_10_New=subset(CPI_Colitis_10_New,subset=nFeature_RNA>900 & nFeature_RNA<1700 & MTpercent<5 & nCount_RNA>2000 & nCount_RNA<4500)
################################################################################ End CPI_Colitis_10_New



################################################################################ Start CPI_Colitis_11_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_11_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_11_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 4500, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_11_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 900, color = 'red') +
  geom_vline(xintercept = 1700, color = 'red')
dev.off()


# CPI_Colitis_11_New=subset(CPI_Colitis_11_New,subset=nFeature_RNA>900 & nFeature_RNA<1700 & MTpercent<5 & nCount_RNA>2000 & nCount_RNA<4500)
################################################################################ End CPI_Colitis_11_New



################################################################################ Start CPI_Colitis_12_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_12_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_12_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 4500, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_12_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 900, color = 'red') +
  geom_vline(xintercept = 1700, color = 'red')
dev.off()


# CPI_Colitis_12_New=subset(CPI_Colitis_12_New,subset=nFeature_RNA>900 & nFeature_RNA<1700 & MTpercent<5 & nCount_RNA>2000 & nCount_RNA<4500)
################################################################################ End CPI_Colitis_12_New



################################################################################ Start CPI_Colitis_13_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_13_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_13_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 4500, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_13_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 900, color = 'red') +
  geom_vline(xintercept = 1700, color = 'red')
dev.off()


# CPI_Colitis_13_New=subset(CPI_Colitis_13_New,subset=nFeature_RNA>900 & nFeature_RNA<1700 & MTpercent<5 & nCount_RNA>2000 & nCount_RNA<4500)
################################################################################ End CPI_Colitis_13_New



################################################################################ Start CPI_Colitis_14_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_14_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_14_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 4500, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_14_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 900, color = 'red') +
  geom_vline(xintercept = 1700, color = 'red')
dev.off()


# CPI_Colitis_14_New=subset(CPI_Colitis_14_New,subset=nFeature_RNA>900 & nFeature_RNA<1700 & MTpercent<5 & nCount_RNA>2000 & nCount_RNA<4500)
################################################################################ End CPI_Colitis_14_New



################################################################################ Start CPI_Colitis_15_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_15_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_15_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 4500, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_15_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 900, color = 'red') +
  geom_vline(xintercept = 1700, color = 'red')
dev.off()


# CPI_Colitis_15_New=subset(CPI_Colitis_15_New,subset=nFeature_RNA>900 & nFeature_RNA<1700 & MTpercent<5 & nCount_RNA>2000 & nCount_RNA<4500)
################################################################################ End CPI_Colitis_15_New



################################################################################ Start CPI_Colitis_16_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_16_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_16_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 4500, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_16_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 900, color = 'red') +
  geom_vline(xintercept = 1700, color = 'red')
dev.off()


# CPI_Colitis_16_New=subset(CPI_Colitis_16_New,subset=nFeature_RNA>900 & nFeature_RNA<1700 & MTpercent<5 & nCount_RNA>2000 & nCount_RNA<4500)
################################################################################ End CPI_Colitis_16_New



################################################################################ Start CPI_Colitis_17_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_17_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_17_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 4500, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_17_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 900, color = 'red') +
  geom_vline(xintercept = 1700, color = 'red')
dev.off()


# CPI_Colitis_17_New=subset(CPI_Colitis_17_New,subset=nFeature_RNA>900 & nFeature_RNA<1700 & MTpercent<5 & nCount_RNA>2000 & nCount_RNA<4500)
################################################################################ End CPI_Colitis_17_New



################################################################################ Start CPI_Colitis_18_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_18_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_18_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 4500, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_18_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 900, color = 'red') +
  geom_vline(xintercept = 1700, color = 'red')
dev.off()


# CPI_Colitis_18_New=subset(CPI_Colitis_18_New,subset=nFeature_RNA>900 & nFeature_RNA<1700 & MTpercent<5 & nCount_RNA>2000 & nCount_RNA<4500)
################################################################################ End CPI_Colitis_18_New



################################################################################ Start CPI_Colitis_19_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_19_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_19_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 5000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_19_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 950, color = 'red') +
  geom_vline(xintercept = 2000, color = 'red')
dev.off()


CPI_Colitis_19_New=subset(CPI_Colitis_19_New,subset=nFeature_RNA>950 & nFeature_RNA<2000 & MTpercent<5 & nCount_RNA>2000 & nCount_RNA<5000)
################################################################################ End CPI_Colitis_19_New



################################################################################ Start CPI_Colitis_20_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_20_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_20_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 5200, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_20_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 950, color = 'red') +
  geom_vline(xintercept = 1800, color = 'red')
dev.off()


CPI_Colitis_20_New=subset(CPI_Colitis_20_New,subset=nFeature_RNA>950 & nFeature_RNA<1800 & MTpercent<5 & nCount_RNA>2000 & nCount_RNA<5200)
################################################################################ End CPI_Colitis_20_New



################################################################################ Start CPI_Colitis_21_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Colitis_21_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_21_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 5000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Colitis_21_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 900, color = 'red') +
  geom_vline(xintercept = 1800, color = 'red')
dev.off()


CPI_Colitis_21_New=subset(CPI_Colitis_21_New,subset=nFeature_RNA>900 & nFeature_RNA<1800 & MTpercent<5 & nCount_RNA>2000 & nCount_RNA<5000)
################################################################################ End CPI_Colitis_21_New



################################################################################ Start CPI_Control_1_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_1_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_1_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 2.5, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 2800, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_1_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 2.5, color = 'red') + 
  geom_vline(xintercept = 350, color = 'red') +
  geom_vline(xintercept = 1300, color = 'red')
dev.off()


CPI_Control_1_New=subset(CPI_Control_1_New,subset=nFeature_RNA>350 & nFeature_RNA<1300 & MTpercent<2.5 & nCount_RNA>600 & nCount_RNA<2800)
################################################################################ End CPI_Control_1_New



################################################################################ Start CPI_Control_2_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_2_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_2_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 3, color = 'red') + 
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_2_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 3, color = 'red') + 
  geom_vline(xintercept = 400, color = 'red') +
  geom_vline(xintercept = 1200, color = 'red')
dev.off()


CPI_Control_2_New=subset(CPI_Control_2_New,subset=nFeature_RNA>400 & nFeature_RNA<1200 & MTpercent<3 & nCount_RNA>700 & nCount_RNA<3000)
################################################################################ End CPI_Control_2_New



################################################################################ Start CPI_Control_3_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_3_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_3_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 3, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_3_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 3, color = 'red') + 
  geom_vline(xintercept = 350, color = 'red') +
  geom_vline(xintercept = 1200, color = 'red')
dev.off()


CPI_Control_3_New=subset(CPI_Control_3_New,subset=nFeature_RNA>350 & nFeature_RNA<1200 & MTpercent<3 & nCount_RNA>600 & nCount_RNA<3000)
################################################################################ End CPI_Control_3_New



################################################################################ Start CPI_Control_4_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_4_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_4_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 3, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 2800, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_4_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 3, color = 'red') + 
  geom_vline(xintercept = 350, color = 'red') +
  geom_vline(xintercept = 1000, color = 'red')
dev.off()


CPI_Control_4_New=subset(CPI_Control_4_New,subset=nFeature_RNA>350 & nFeature_RNA<1000 & MTpercent<3 & nCount_RNA>600 & nCount_RNA<2800)
################################################################################ End CPI_Control_4_New



################################################################################ Start CPI_Control_5_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_5_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_5_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 2.5, color = 'red') + 
  geom_vline(xintercept = 550, color = 'red') +
  geom_vline(xintercept = 2900, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_5_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 2.5, color = 'red') + 
  geom_vline(xintercept = 350, color = 'red') +
  geom_vline(xintercept = 1200, color = 'red')
dev.off()


CPI_Control_5_New=subset(CPI_Control_5_New,subset=nFeature_RNA>350 & nFeature_RNA<1200 & MTpercent<2.5 & nCount_RNA>550 & nCount_RNA<2900)
################################################################################ End CPI_Control_5_New



################################################################################ Start CPI_Control_6_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_6_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_6_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_6_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 1000, color = 'red')
dev.off()


CPI_Control_6_New=subset(CPI_Control_6_New,subset=nFeature_RNA>600 & nFeature_RNA<1000 & MTpercent<5 & nCount_RNA>1200 & nCount_RNA<3000)
################################################################################ End CPI_Control_6_New



################################################################################ Start CPI_Control_7_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_7_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_7_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_7_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 1400, color = 'red')
dev.off()


CPI_Control_7_New=subset(CPI_Control_7_New,subset=nFeature_RNA>600 & nFeature_RNA<1400 & MTpercent<4 & nCount_RNA>1200 & nCount_RNA<3000)
################################################################################ End CPI_Control_7_New



################################################################################ Start CPI_Control_8_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_8_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_8_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_8_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 1400, color = 'red')
dev.off()


CPI_Control_8_New=subset(CPI_Control_8_New,subset=nFeature_RNA>600 & nFeature_RNA<1400 & MTpercent<4 & nCount_RNA>1200 & nCount_RNA<3000)
################################################################################ End CPI_Control_8_New



################################################################################ Start CPI_Control_9_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_9_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_9_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_9_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 1400, color = 'red')
dev.off()


CPI_Control_9_New=subset(CPI_Control_9_New,subset=nFeature_RNA>600 & nFeature_RNA<1400 & MTpercent<4 & nCount_RNA>1200 & nCount_RNA<3000)
################################################################################ End CPI_Control_9_New



################################################################################ Start CPI_Control_10_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_10_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_10_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_10_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 1400, color = 'red')
dev.off()


CPI_Control_10_New=subset(CPI_Control_10_New,subset=nFeature_RNA>600 & nFeature_RNA<1400 & MTpercent<4 & nCount_RNA>1200 & nCount_RNA<3000)
################################################################################ End CPI_Control_10_New



################################################################################ Start CPI_Control_11_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_11_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_11_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_11_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 1400, color = 'red')
dev.off()


CPI_Control_11_New=subset(CPI_Control_11_New,subset=nFeature_RNA>600 & nFeature_RNA<1400 & MTpercent<4 & nCount_RNA>1200 & nCount_RNA<3000)
################################################################################ End CPI_Control_11_New




################################################################################ Start CPI_Control_12_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_12_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_12_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_12_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 1400, color = 'red')
dev.off()


CPI_Control_12_New=subset(CPI_Control_12_New,subset=nFeature_RNA>600 & nFeature_RNA<1400 & MTpercent<4 & nCount_RNA>1200 & nCount_RNA<3000)
################################################################################ End CPI_Control_12_New



################################################################################ Start CPI_Control_13_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_13_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_13_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_13_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 1400, color = 'red')
dev.off()


CPI_Control_13_New=subset(CPI_Control_13_New,subset=nFeature_RNA>600 & nFeature_RNA<1400 & MTpercent<4 & nCount_RNA>1200 & nCount_RNA<3000)
################################################################################ End CPI_Control_13_New



################################################################################ Start CPI_Control_14_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_14_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_14_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_14_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 1400, color = 'red')
dev.off()


CPI_Control_14_New=subset(CPI_Control_14_New,subset=nFeature_RNA>600 & nFeature_RNA<1400 & MTpercent<4 & nCount_RNA>1200 & nCount_RNA<3000)
################################################################################ End CPI_Control_14_New



################################################################################ Start CPI_Control_15_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_15_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_15_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_15_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 1400, color = 'red')
dev.off()


CPI_Control_15_New=subset(CPI_Control_15_New,subset=nFeature_RNA>600 & nFeature_RNA<1400 & MTpercent<4 & nCount_RNA>1200 & nCount_RNA<3000)
################################################################################ End CPI_Control_15_New



################################################################################ Start CPI_Control_16_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/CPI_Control_16_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_16_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = CPI_Control_16_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 1400, color = 'red')
dev.off()


CPI_Control_16_New=subset(CPI_Control_16_New,subset=nFeature_RNA>600 & nFeature_RNA<1400 & MTpercent<4 & nCount_RNA>1200 & nCount_RNA<3000)
################################################################################ End CPI_Control_16_New



################################################################################ Start Healthy_Control_1_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/Healthy_Control_1_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_1_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 5000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_1_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6, color = 'red') + 
  geom_vline(xintercept = 800, color = 'red') +
  geom_vline(xintercept = 1600, color = 'red')
dev.off()


Healthy_Control_1_New=subset(Healthy_Control_1_New,subset=nFeature_RNA>800 & nFeature_RNA<1600 & MTpercent<6 & nCount_RNA>2000 & nCount_RNA<5000)
################################################################################ End Healthy_Control_1_New



################################################################################ Start Healthy_Control_2_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/Healthy_Control_2_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_2_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6.5, color = 'red') + 
  geom_vline(xintercept = 1300, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_2_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6.5, color = 'red') + 
  geom_vline(xintercept = 550, color = 'red') +
  geom_vline(xintercept = 1000, color = 'red')
dev.off()


Healthy_Control_2_New=subset(Healthy_Control_2_New,subset=nFeature_RNA>550 & nFeature_RNA<1000 & MTpercent<6.5 & nCount_RNA>1300 & nCount_RNA<3000)
################################################################################ End Healthy_Control_2_New



################################################################################ Start Healthy_Control_3_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/Healthy_Control_3_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_3_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6.5, color = 'red') + 
  geom_vline(xintercept = 1300, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_3_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6.5, color = 'red') + 
  geom_vline(xintercept = 550, color = 'red') +
  geom_vline(xintercept = 1000, color = 'red')
dev.off()


# Healthy_Control_3_New=subset(Healthy_Control_3_New,subset=nFeature_RNA>550 & nFeature_RNA<1000 & MTpercent<6.5 & nCount_RNA>1300 & nCount_RNA<3000)
################################################################################ End Healthy_Control_3_New



################################################################################ Start Healthy_Control_4_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/Healthy_Control_4_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_4_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6.5, color = 'red') + 
  geom_vline(xintercept = 1300, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_4_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6.5, color = 'red') + 
  geom_vline(xintercept = 550, color = 'red') +
  geom_vline(xintercept = 1000, color = 'red')
dev.off()


# Healthy_Control_4_New=subset(Healthy_Control_4_New,subset=nFeature_RNA>550 & nFeature_RNA<1000 & MTpercent<6.5 & nCount_RNA>1300 & nCount_RNA<3000)
################################################################################ End Healthy_Control_4_New



################################################################################ Start Healthy_Control_5_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/Healthy_Control_5_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_5_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6, color = 'red') + 
  geom_vline(xintercept = 1500, color = 'red') +
  geom_vline(xintercept = 3300, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_5_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6, color = 'red') + 
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 1400, color = 'red')
dev.off()


Healthy_Control_5_New=subset(Healthy_Control_5_New,subset=nFeature_RNA>700 & nFeature_RNA<1400 & MTpercent<6 & nCount_RNA>1500 & nCount_RNA<3300)
################################################################################ End Healthy_Control_5_New



################################################################################ Start Healthy_Control_6_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/Healthy_Control_6_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_6_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6, color = 'red') + 
  geom_vline(xintercept = 1500, color = 'red') +
  geom_vline(xintercept = 3300, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_6_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6, color = 'red') + 
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 1400, color = 'red')
dev.off()


# Healthy_Control_6_New=subset(Healthy_Control_6_New,subset=nFeature_RNA>700 & nFeature_RNA<1400 & MTpercent<6 & nCount_RNA>1500 & nCount_RNA<3300)
################################################################################ End Healthy_Control_6_New



################################################################################ Start Healthy_Control_7_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/Healthy_Control_7_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_7_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6.5, color = 'red') + 
  geom_vline(xintercept = 1800, color = 'red') +
  geom_vline(xintercept = 5000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_7_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6.5, color = 'red') + 
  geom_vline(xintercept = 800, color = 'red') +
  geom_vline(xintercept = 1500, color = 'red')
dev.off()


Healthy_Control_7_New=subset(Healthy_Control_7_New,subset=nFeature_RNA>800 & nFeature_RNA<1500 & MTpercent<6.5 & nCount_RNA>1800 & nCount_RNA<5000)
################################################################################ End Healthy_Control_7_New



################################################################################ Start Healthy_Control_8_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/Healthy_Control_8_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_8_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6.5, color = 'red') + 
  geom_vline(xintercept = 1900, color = 'red') +
  geom_vline(xintercept = 5000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_8_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6.5, color = 'red') + 
  geom_vline(xintercept = 900, color = 'red') +
  geom_vline(xintercept = 1600, color = 'red')
dev.off()


Healthy_Control_8_New=subset(Healthy_Control_8_New,subset=nFeature_RNA>900 & nFeature_RNA<1600 & MTpercent<6.5 & nCount_RNA>1900 & nCount_RNA<5000)
################################################################################ End Healthy_Control_8_New



################################################################################ Start Healthy_Control_9_New
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3/Quality Control/Healthy_Control_9_New")


png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_9_New@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6, color = 'red') + 
  geom_vline(xintercept = 1850, color = 'red') +
  geom_vline(xintercept = 6000, color = 'red')
dev.off()


png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Healthy_Control_9_New@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 6, color = 'red') + 
  geom_vline(xintercept = 800, color = 'red') +
  geom_vline(xintercept = 1800, color = 'red')
dev.off()


Healthy_Control_9_New=subset(Healthy_Control_9_New,subset=nFeature_RNA>800 & nFeature_RNA<1800 & MTpercent<6 & nCount_RNA>1850 & nCount_RNA<6000)
################################################################################ End Healthy_Control_9_New