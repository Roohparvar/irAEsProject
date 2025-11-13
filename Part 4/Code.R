################################################################################ Start libraries
library(Seurat)
library(SeuratObject)
library(ggplot2)
################################################################################ End libraries



################################################################################ Start reading the expression and metadata CSVs, and splitting the data by each patient.
# Expression matrix
expr <- read.csv("C:/Esmaeil/irAEsProject/Backup/Part 4/0_Data/Python Data Preparation/expression.csv", row.names = 1)

# Metadata file
obs <- read.csv("C:/Esmaeil/irAEsProject/Backup/Part 4/0_Data/Python Data Preparation/metadata_obs.csv", row.names = 1)



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

# Saved and stored: HS7 (3172 cells, disease: CPI_colitis)
# Saved and stored: HS8 (2851 cells, disease: CPI_colitis)
# Saved and stored: HS9 (2366 cells, disease: CPI_colitis)
# Saved and stored: HS10 (2731 cells, disease: CPI_colitis)
# Saved and stored: HS11 (4848 cells, disease: CPI_colitis)
# Saved and stored: HS12 (3589 cells, disease: CPI_colitis)
# Saved and stored: HS13 (5142 cells, disease: CPI_colitis)
# Saved and stored: HS14 (5293 cells, disease: CPI_colitis)
# Saved and stored: HS15 (4936 cells, disease: CPI_colitis)
################################################################################ End reading the expression and metadata CSVs, and splitting the data by each patient.

# QC

################################################################################ Start HS1(HC)
HS1 <- sample_list[["HS1"]]
CountMatrix <- t(HS1$expression)  


Jun_srobj_1 <- CreateSeuratObject(counts = CountMatrix, project = "Jun_HS1_HC", min.cells = 3, min.features = 200)

# grep("^MT", rownames(Jun_srobj_1), value = TRUE)
Jun_srobj_1[["MTpercent"]] <- PercentageFeatureSet(Jun_srobj_1, pattern = "^MT\\.")

setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4/QC/HS1(HC)")

png(filename = "1.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_1@meta.data, aes(x = nCount_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 15, color = 'red') +
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 30000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_1@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 15, color = 'red') +
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


Jun_srobj_1 <- subset(Jun_srobj_1, subset = nFeature_RNA > 600 & nFeature_RNA < 3000 & MTpercent < 15 & nCount_RNA > 1000 & nCount_RNA < 30000)
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
  geom_hline(yintercept = 18, color = 'red') +
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 30000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_2@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 18, color = 'red') +
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


Jun_srobj_2 <- subset(Jun_srobj_2, subset = nFeature_RNA > 600 & nFeature_RNA < 3000 & MTpercent < 18 & nCount_RNA > 1000 & nCount_RNA < 30000)
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
  geom_hline(yintercept = 18, color = 'red') +
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 30000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_3@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 18, color = 'red') +
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


Jun_srobj_3 <- subset(Jun_srobj_3, subset = nFeature_RNA > 700 & nFeature_RNA < 3000 & MTpercent < 18 & nCount_RNA > 1200 & nCount_RNA < 30000)
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
  geom_hline(yintercept = 15, color = 'red') +
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 30000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_4@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 15, color = 'red') +
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


Jun_srobj_4 <- subset(Jun_srobj_4, subset = nFeature_RNA > 700 & nFeature_RNA < 3000 & MTpercent < 15 & nCount_RNA > 1200 & nCount_RNA < 30000)
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
  geom_hline(yintercept = 18, color = 'red') + 
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 7900, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_5@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 18, color = 'red') + 
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 2500, color = 'red')
dev.off()


Jun_srobj_5 <- subset(Jun_srobj_5, subset = nFeature_RNA > 700 & nFeature_RNA < 2500 & MTpercent < 18 & nCount_RNA > 1000 & nCount_RNA < 7900)
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
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 10000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_7@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


Jun_srobj_7 <- subset(Jun_srobj_7, subset = nFeature_RNA > 700 & nFeature_RNA < 3000 & MTpercent < 20 & nCount_RNA > 1000 & nCount_RNA < 10000)
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
  geom_hline(yintercept = 17, color = 'red') + 
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 35000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_8@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 17, color = 'red') + 
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 3300, color = 'red')
dev.off()


Jun_srobj_8 <- subset(Jun_srobj_8, subset = nFeature_RNA > 700 & nFeature_RNA < 3300 & MTpercent < 17 & nCount_RNA > 1200 & nCount_RNA < 35000)
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
  geom_hline(yintercept = 18, color = 'red') + 
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 14000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_9@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 18, color = 'red') + 
  geom_vline(xintercept = 650, color = 'red') +
  geom_vline(xintercept = 3700, color = 'red')
dev.off()


Jun_srobj_9 <- subset(Jun_srobj_9, subset = nFeature_RNA > 650 & nFeature_RNA < 3700 & MTpercent < 18 & nCount_RNA > 1200 & nCount_RNA < 14000)
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
  geom_hline(yintercept = 18, color = 'red') + 
  geom_vline(xintercept = 1400, color = 'red') +
  geom_vline(xintercept = 15000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_10@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 18, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


Jun_srobj_10 <- subset(Jun_srobj_10, subset = nFeature_RNA > 600 & nFeature_RNA < 3000 & MTpercent < 18 & nCount_RNA > 1400 & nCount_RNA < 15000)
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
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 18000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_11@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 3200, color = 'red')
dev.off()


Jun_srobj_11 <- subset(Jun_srobj_11, subset = nFeature_RNA > 600 & nFeature_RNA < 3200 & MTpercent < 20 & nCount_RNA > 1000 & nCount_RNA < 18000)
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
  geom_hline(yintercept = 18, color = 'red') + 
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 10000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_12@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 18, color = 'red') + 
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


Jun_srobj_12 <- subset(Jun_srobj_12, subset = nFeature_RNA > 700 & nFeature_RNA < 3000 & MTpercent < 18 & nCount_RNA > 1000 & nCount_RNA < 10000)
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
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 17000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_13@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 3200, color = 'red')
dev.off()


Jun_srobj_13 <- subset(Jun_srobj_13, subset = nFeature_RNA > 600 & nFeature_RNA < 3200 & MTpercent < 20 & nCount_RNA > 1000 & nCount_RNA < 17000)
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
  geom_hline(yintercept = 17, color = 'red') + 
  geom_vline(xintercept = 1100, color = 'red') +
  geom_vline(xintercept = 42000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_14@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 17, color = 'red') + 
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 8000, color = 'red')
dev.off()


Jun_srobj_14 <- subset(Jun_srobj_14, subset = nFeature_RNA > 700 & nFeature_RNA < 8000 & MTpercent < 17 & nCount_RNA > 1100 & nCount_RNA < 42000)
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
  geom_hline(yintercept = 16, color = 'red') + 
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 40000, color = 'red')
dev.off()

png(filename = "2.png", width = 10000, height = 4000, units = "px", res = 600)
ggplot(data = Jun_srobj_15@meta.data, aes(x = nFeature_RNA, y = MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nFeature_RNA', y = 'Percent Mito') +
  scale_x_log10() +
  geom_hline(yintercept = 16, color = 'red') + 
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 5400, color = 'red')
dev.off()


Jun_srobj_15 <- subset(Jun_srobj_15, subset = nFeature_RNA > 700 & nFeature_RNA < 5400 & MTpercent < 16 & nCount_RNA > 1000 & nCount_RNA < 40000)
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


setwd("C:/Esmaeil/irAEsProject/Backup/Part 4/1_The Seurat object obtained after RunPCA and before IntegrateLayers")
# 1
# saveRDS(file = "merged_obj",merged_obj)
# The Seurat object obtained after RunPCA and before IntegrateLayers


merged_obj <- IntegrateLayers(object = merged_obj,
                              method = CCAIntegration,
                              orig.reduction = "pca", 
                              new.reduction = "integrated.cca",
                              verbose = FALSE)

# 2
setwd("C:/Esmaeil/irAEsProject/Backup/Part 4/2_The Seurat object obtained after IntegrateLayers")
# saveRDS(file = "merged_obj",merged_obj)
# The Seurat object obtained after IntegrateLayers

merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])
################################################################################ End integration

merged_obj1 = merged_obj

################################################################################ Start UMAP
merged_obj1=FindNeighbors(merged_obj1,dims = 1:30,reduction = "integrated.cca")
merged_obj1=FindClusters(merged_obj1,resolution = 0.1)


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 4")

merged_obj1=RunUMAP(merged_obj1,dims = 1:30, reduction = "integrated.cca")
png(filename = "1_First UMAP.png", width = 4000,height = 3000, units ="px",res = 600 )
DimPlot(merged_obj1,label = TRUE)
dev.off()



png(filename = "1_First Dot Plot.png", width = 20000, height = 4000, units = "px", res = 600)


DotPlot(merged_obj1, features = c(
  "CD79A", "MS4A1", "CD19",
  "CD3D", "IL7R", "CD4", "CD40LG",
  "CD8A", "CD8B", "GZMA",
  "MKI67", "TYMS", "STMN1",
  "PLVAP", "RAMP2", "EGFL7", "PECAM1", "CLDN5", "CRYAB", "S100B",
  "NRXN1", "PLP1", "FABP1", "LRT8", "LGALS4", "PIGR", "EPCAM",
  "TRDC", "TRGC2", "TRGC1",
  "KRT86", "KRT81", "RORC",
  "TPSAB1", "CPA3", "GATA2",
  "COL3A1", "COL1A2", "PDGFRA",
  "LYZ", "C1QA", "C1QB", "C1QC",
  "GNLY", "KLRF1", "FCGR3A",
  "MZB1", "DERL3", "TNFRSF17", "JCHAIN"
)) +
  ggplot2::theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10)
  ) +
  ggplot2::labs(
    title = "Cluster Marker Gene Expression (Dot Plot)",
    x = "Genes",
    y = "Clusters"
  ) +
  ggplot2::scale_size(range = c(1, 8)) 

dev.off()



# Remove clusters 1, 3 ,4, 6, 7, 8, 9, 10, and 11
merged_obj1 <- subset(merged_obj1, subset = seurat_clusters %in% c(1, 3, 4, 6, 7, 8, 9, 10, 11), invert = TRUE)

################################################################################ End UMAP

merged_obj2 = merged_obj1

merged_obj2$orig.ident <- gsub("^Jun_HS([0-9]+)_HC$", "Healthy HS\\1", merged_obj2$orig.ident)
merged_obj2$orig.ident <- gsub("^Jun_HS([0-9]+)_UC$", "UC_Inflamed HS\\1", merged_obj2$orig.ident)
merged_obj2$orig.ident <- gsub("^Jun_HS([0-9]+)_CPI_colitis$", "CPI_Colitis HS\\1", merged_obj2$orig.ident)

merged_obj2 <- RenameCells(merged_obj2, new.names = paste0("Jun/", merged_obj2$orig.ident, "/", colnames(merged_obj2)))

cell_counts <- as.data.frame(table(merged_obj2$orig.ident))
colnames(cell_counts) <- c("Sample", "Number_of_Cells")
print(cell_counts, row.names = FALSE)

################################################################################ Start Extracting and saving Seurat objects for each sample
samples <- unique(merged_obj2$orig.ident)

start_index <- 84
num_samples <- 14 


for (i in 1:num_samples) {
  sample_name <- samples[i] 
  seurat_obj <- subset(merged_obj2, subset = orig.ident == sample_name)
  assign(paste0("srobj_", start_index + i - 1), seurat_obj)  # srobj_84 ... srobj_97
}


# 3
# 3_The Seurat objects per sample
setwd("C:/Esmaeil/irAEsProject/Backup/Part 4/3_The Seurat objects per sample")

seurat_objs <- mget(paste0("srobj_", start_index:(start_index + num_samples - 1)))
for (i in seq_along(seurat_objs)) {
  saveRDS(seurat_objs[[i]], file = paste0("srobj_", start_index + i - 1, ".rds"))
}

################################################################################ End Extracting and saving Seurat objects for each sample