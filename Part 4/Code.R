################################################################################ Start libraries
library(Seurat)
library(SeuratObject)
library(ggplot2)
################################################################################ End libraries



################################################################################ Start reading the expression and metadata CSVs, and splitting the data by each patient.
# Expression matrix
expr <- read.csv("C:/Esmaeil/irAEsProject/Backup/Part 4/Data/1_Python Data Preparation/expression.csv", row.names = 1)

# Metadata file
obs <- read.csv("C:/Esmaeil/irAEsProject/Backup/Part 4/Data/1_Python Data Preparation/metadata_obs.csv", row.names = 1)



# Look at the column names to find the group info
colnames(obs)

# Unique samples (based on our data)
samples <- c("HS1","HS2","HS3","HS4","HS5","HS7","HS8","HS9", "HS10","HS11","HS12","HS13","HS14","HS15")

# Check if rownames match
all(rownames(expr) == rownames(obs))

# Convert to character if needed
obs$CoLabs_patient <- as.character(obs$CoLabs_patient)



# Keep results in memory
sample_list <- list()

for (s in samples) {
  cells <- rownames(obs)[obs$CoLabs_patient == s]
  
  expr_subset <- expr[cells, , drop = FALSE]
  obs_subset  <- obs[cells, , drop = FALSE]
  
  sample_data <- list(expression = expr_subset, metadata = obs_subset)
  
  sample_list[[s]] <- sample_data
  
  cat(s, "(", nrow(expr_subset), "cells, disease:",
      unique(obs_subset$disease), ")\n")
}

# Saved and stored: HS1 ( 3542 cells, disease: HC )
# Saved and stored: HS2 ( 3286 cells, disease: HC )
# Saved and stored: HS3 ( 3976 cells, disease: HC )

# Saved and stored: HS4 ( 2853 cells, disease: UC )
# Saved and stored: HS5 ( 2260 cells, disease: UC )

# Saved and stored: HS7 ( 3172 cells, disease: CPI_colitis )
# Saved and stored: HS8 ( 2851 cells, disease: CPI_colitis )
# Saved and stored: HS9 ( 2366 cells, disease: CPI_colitis )
# Saved and stored: HS10 ( 2731 cells, disease: CPI_colitis )
# Saved and stored: HS11 ( 4848 cells, disease: CPI_colitis )
# Saved and stored: HS12 ( 3589 cells, disease: CPI_colitis )
# Saved and stored: HS13 ( 5142 cells, disease: CPI_colitis )
# Saved and stored: HS14 ( 5293 cells, disease: CPI_colitis )
# Saved and stored: HS15 ( 4936 cells, disease: CPI_colitis )
################################################################################ End reading the expression and metadata CSVs, and splitting the data by each patient.

# QC

################################################################################ Start HS1(HC)
HS1 <- sample_list[["HS1"]]
CountMatrix <- t(HS1$expression)  


Jun_srobj_1 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS1_HC", min.cells = 3, min.features = 200)
Jun_srobj_1[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_1, pattern = "^MT\\.")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS1(HC)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_1@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 3500, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_1@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 650, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


Jun_srobj_1 <- subset(Jun_srobj_1, subset = nFeature_RNA > 650 & nFeature_RNA < 3000 & MTpercent < 5 & nCount_RNA > 1000 & nCount_RNA < 3500)
################################################################################ End HS1(HC)


################################################################################ Start HS2(HC)
HS2 <- sample_list[["HS2"]]
CountMatrix <- t(HS2$expression)  


Jun_srobj_2 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS2_HC", min.cells = 3, min.features = 200)
Jun_srobj_2[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_2, pattern = "^MT\\.")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS2(HC)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_2@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 3500, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_2@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 650, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


Jun_srobj_2 <- subset(Jun_srobj_2, subset = nFeature_RNA > 650 & nFeature_RNA < 3000 & MTpercent < 5 & nCount_RNA > 1000 & nCount_RNA < 3500)
################################################################################ End HS2(HC)



################################################################################ Start HS3(HC)
HS3 <- sample_list[["HS3"]]
CountMatrix <- t(HS3$expression)  


Jun_srobj_3 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS3_HC", min.cells = 3, min.features = 200)
Jun_srobj_3[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_3, pattern = "^MT\\.")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS3(HC)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_3@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 3500, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_3@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 650, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


Jun_srobj_3 <- subset(Jun_srobj_3, subset = nFeature_RNA > 650 & nFeature_RNA < 3000 & MTpercent < 5 & nCount_RNA > 1000 & nCount_RNA < 3500)
################################################################################ End HS3(HC)


################################################################################ Start HS4(UC)
HS4 <- sample_list[["HS4"]]
CountMatrix <- t(HS4$expression)  


Jun_srobj_4 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS4_UC", min.cells = 3, min.features = 200)
Jun_srobj_4[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_4, pattern = "^MT\\.")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS4(UC)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_4@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_4@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 650, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


Jun_srobj_4 <- subset(Jun_srobj_4, subset = nFeature_RNA > 650 & nFeature_RNA < 3000 & MTpercent < 5 & nCount_RNA > 1000 & nCount_RNA < 3000)
################################################################################ End HS4(UC)


################################################################################ Start HS5(UC)
HS5 <- sample_list[["HS5"]]
CountMatrix <- t(HS5$expression)  


Jun_srobj_5 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS5_UC", min.cells = 3, min.features = 200)
Jun_srobj_5[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_5, pattern = "^MT\\.")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS5(UC)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_5@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 1500, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_5@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 2800, color = 'red')
dev.off()


Jun_srobj_5 <- subset(Jun_srobj_5, subset = nFeature_RNA > 600 & nFeature_RNA < 2800 & MTpercent < 5 & nCount_RNA > 1500 & nCount_RNA < 3000)
################################################################################ End HS5(UC)



################################################################################ Start HS7(CPI_colitis)
HS7 <- sample_list[["HS7"]]
CountMatrix <- t(HS7$expression)  


Jun_srobj_7 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS7_CPI_colitis", min.cells = 3, min.features = 200)
Jun_srobj_7[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_7, pattern = "^MT\\.")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS7(CPI_colitis)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_7@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4.5, color = 'red') +
  geom_vline(xintercept = 1500, color = 'red') +
  geom_vline(xintercept = 3500, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_7@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4.5, color = 'red') +
  geom_vline(xintercept = 500, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


Jun_srobj_7 <- subset(Jun_srobj_7, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & MTpercent < 4.5 & nCount_RNA > 1500 & nCount_RNA < 3500)
################################################################################ End HS7(CPI_colitis)



################################################################################ Start HS8(CPI_colitis)
HS8 <- sample_list[["HS8"]]
CountMatrix <- t(HS8$expression)  


Jun_srobj_8 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS8_CPI_colitis", min.cells = 3, min.features = 200)
Jun_srobj_8[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_8, pattern = "^MT\\.")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS8(CPI_colitis)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_8@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 3500, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_8@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 500, color = 'red') +
  geom_vline(xintercept = 3100, color = 'red')
dev.off()


Jun_srobj_8 <- subset(Jun_srobj_8, subset = nFeature_RNA > 500 & nFeature_RNA < 3100 & MTpercent < 5 & nCount_RNA > 1000 & nCount_RNA < 3500)
################################################################################ End HS8(CPI_colitis)



################################################################################ Start HS9(CPI_colitis)
HS9 <- sample_list[["HS9"]]
CountMatrix <- t(HS9$expression)  


Jun_srobj_9 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS9_CPI_colitis", min.cells = 3, min.features = 200)
Jun_srobj_9[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_9, pattern = "^MT\\.")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS9(CPI_colitis)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_9@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') +
  geom_vline(xintercept = 1500, color = 'red') +
  geom_vline(xintercept = 3600, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_9@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4, color = 'red') +
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 4000, color = 'red')
dev.off()


Jun_srobj_9 <- subset(Jun_srobj_9, subset = nFeature_RNA > 600 & nFeature_RNA < 4000 & MTpercent < 4 & nCount_RNA > 1500 & nCount_RNA < 3600)
################################################################################ End HS9(CPI_colitis)


################################################################################ Start HS10(CPI_colitis)
HS10 <- sample_list[["HS10"]]
CountMatrix <- t(HS10$expression)  


Jun_srobj_10 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS10_CPI_colitis", min.cells = 3, min.features = 200)
Jun_srobj_10[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_10, pattern = "^MT\\.")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS10(CPI_colitis)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_10@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 3700, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_10@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 500, color = 'red') +
  geom_vline(xintercept = 4000, color = 'red')
dev.off()


Jun_srobj_10 <- subset(Jun_srobj_10, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & MTpercent < 5 & nCount_RNA > 1200 & nCount_RNA < 3700)
################################################################################ End HS10(CPI_colitis)



################################################################################ Start HS11(CPI_colitis)
HS11 <- sample_list[["HS11"]]
CountMatrix <- t(HS11$expression)  


Jun_srobj_11 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS11_CPI_colitis", min.cells = 3, min.features = 200)
Jun_srobj_11[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_11, pattern = "^MT\\.")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS11(CPI_colitis)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_11@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4.5, color = 'red') +
  geom_vline(xintercept = 1300, color = 'red') +
  geom_vline(xintercept = 3700, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_11@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4.5, color = 'red') +
  geom_vline(xintercept = 500, color = 'red') +
  geom_vline(xintercept = 4000, color = 'red')
dev.off()


Jun_srobj_11 <- subset(Jun_srobj_11, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & MTpercent < 4.5 & nCount_RNA > 1300 & nCount_RNA < 3700)
################################################################################ End HS11(CPI_colitis)



################################################################################ Start HS12(CPI_colitis)
HS12 <- sample_list[["HS12"]]
CountMatrix <- t(HS12$expression)  


Jun_srobj_12 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS12_CPI_colitis", min.cells = 3, min.features = 200)
Jun_srobj_12[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_12, pattern = "^MT\\.")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS12(CPI_colitis)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_12@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4.5, color = 'red') +
  geom_vline(xintercept = 1300, color = 'red') +
  geom_vline(xintercept = 3700, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_12@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 4.5, color = 'red') +
  geom_vline(xintercept = 500, color = 'red') +
  geom_vline(xintercept = 4000, color = 'red')
dev.off()


Jun_srobj_12 <- subset(Jun_srobj_12, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & MTpercent < 4.5 & nCount_RNA > 1300 & nCount_RNA < 3700)
################################################################################ End HS12(CPI_colitis)



################################################################################ Start HS13(CPI_colitis)
HS13 <- sample_list[["HS13"]]
CountMatrix <- t(HS13$expression)  


Jun_srobj_13 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS13_CPI_colitis", min.cells = 3, min.features = 200)
Jun_srobj_13[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_13, pattern = "^MT\\.")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS13(CPI_colitis)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_13@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 1300, color = 'red') +
  geom_vline(xintercept = 3700, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_13@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 500, color = 'red') +
  geom_vline(xintercept = 4000, color = 'red')
dev.off()


Jun_srobj_13 <- subset(Jun_srobj_13, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & MTpercent < 5 & nCount_RNA > 1300 & nCount_RNA < 3700)
################################################################################ End HS13(CPI_colitis)



################################################################################ Start HS14(CPI_colitis)
HS14 <- sample_list[["HS14"]]
CountMatrix <- t(HS14$expression)  


Jun_srobj_14 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS14_CPI_colitis", min.cells = 3, min.features = 200)
Jun_srobj_14[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_14, pattern = "^MT\\.")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS14(CPI_colitis)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_14@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 1100, color = 'red') +
  geom_vline(xintercept = 4100, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_14@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 500, color = 'red') +
  geom_vline(xintercept = 8500, color = 'red')
dev.off()


Jun_srobj_14 <- subset(Jun_srobj_14, subset = nFeature_RNA > 500 & nFeature_RNA < 8500 & MTpercent < 5 & nCount_RNA > 1100 & nCount_RNA < 4100)
################################################################################ End HS14(CPI_colitis)



################################################################################ Start HS15(CPI_colitis)
HS15 <- sample_list[["HS15"]]
CountMatrix <- t(HS15$expression)  


Jun_srobj_15 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS15_CPI_colitis", min.cells = 3, min.features = 200)
Jun_srobj_15[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_15, pattern = "^MT\\.")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS15(CPI_colitis)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_15@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 4100, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_15@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 5, color = 'red') +
  geom_vline(xintercept = 500, color = 'red') +
  geom_vline(xintercept = 8000, color = 'red')
dev.off()


Jun_srobj_15 <- subset(Jun_srobj_15, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & MTpercent < 5 & nCount_RNA > 1000 & nCount_RNA < 4100)
################################################################################ End HS15(CPI_colitis)



################################################################################ Start integration and UMAP
merged_obj <- merge(Jun_srobj_1, y = list(Jun_srobj_2, Jun_srobj_3, Jun_srobj_4, 
                                          Jun_srobj_5, Jun_srobj_7, Jun_srobj_8, 
                                          Jun_srobj_9, Jun_srobj_10,Jun_srobj_11, 
                                          Jun_srobj_12, Jun_srobj_13,Jun_srobj_14,
                                          Jun_srobj_15))


merged_obj=NormalizeData(merged_obj,normalization.method = "LogNormalize",scale.factor = 10000)

merged_obj=FindVariableFeatures(merged_obj,selection.method = "vst",nfeatures = 2000)

merged_obj = ScaleData(merged_obj,features = rownames(merged_obj))

merged_obj = RunPCA(merged_obj)


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4")
# 1
# saveRDS(file = "merged_obj",merged_obj)
# The Seurat object obtained after RunPCA and before IntegrateLayers


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