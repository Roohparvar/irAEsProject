################################################################################ Start library
library(Seurat)
library(SeuratObject)
library(ggplot2)

install.packages("BiocManager")
BiocManager::install("rhdf5")
library(rhdf5)
################################################################################ End library



################################################################################ Start Control Healthy MC_1 #1
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - Healthy/GSM6250065_MC_1_Colon_555_CD45pos_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - Healthy/GSM6250065_MC_1_Colon_555_CD45pos_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_1 <- CreateSeuratObject(counts = raw_data, project = "Control Healthy MC_1",min.cells=3,min.features=200)
seurat_obj_1[["MTpercent"]]=PercentageFeatureSet(seurat_obj_1,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_1@meta.data, aes(x = seurat_obj_1$nCount_RNA, y = seurat_obj_1$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 60, color = 'red') + 
  geom_vline(xintercept = 400, color = 'red') +
  geom_vline(xintercept = 20000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_1@meta.data, aes(x=seurat_obj_1$nFeature_RNA, y=seurat_obj_1$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=60, color='red') + 
  geom_vline(xintercept=195, color='red') +
  geom_vline(xintercept=2000, color='red') 
dev.off()


seurat_obj_1=subset(seurat_obj_1,subset=nFeature_RNA>195 & nFeature_RNA<2000 & MTpercent<60 & nCount_RNA>400 & nCount_RNA<20000)
################################################################################ End Control Healthy MC_1 #1



################################################################################ Start Control Healthy MC_2_A #2
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - Healthy/GSM6250066_MC_2_Colon_576_CD45sort_5p_GEX_A")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - Healthy/GSM6250066_MC_2_Colon_576_CD45sort_5p_GEX_A/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_2 <- CreateSeuratObject(counts = raw_data, project = "Control Healthy MC_2_A",min.cells=3,min.features=200)
seurat_obj_2[["MTpercent"]]=PercentageFeatureSet(seurat_obj_2,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_2@meta.data, aes(x = seurat_obj_2$nCount_RNA, y = seurat_obj_2$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 30, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 15000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_2@meta.data, aes(x=seurat_obj_2$nFeature_RNA, y=seurat_obj_2$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=30, color='red') + 
  geom_vline(xintercept=300, color='red') +
  geom_vline(xintercept=2700, color='red') 
dev.off()


seurat_obj_2=subset(seurat_obj_2,subset=nFeature_RNA>300 & nFeature_RNA<2700 & MTpercent<30 & nCount_RNA>600 & nCount_RNA<15000)
################################################################################ End Control Healthy MC_2_A #2



################################################################################ Start Control Healthy MC_2_B #3
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - Healthy/GSM6250067_MC_2_Colon_576_CD45sort_5p_GEX_B")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - Healthy/GSM6250067_MC_2_Colon_576_CD45sort_5p_GEX_B/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_3 <- CreateSeuratObject(counts = raw_data, project = "Control Healthy MC_2_B",min.cells=3,min.features=200)
seurat_obj_3[["MTpercent"]]=PercentageFeatureSet(seurat_obj_3,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_3@meta.data, aes(x = seurat_obj_3$nCount_RNA, y = seurat_obj_3$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 30, color = 'red') + 
  geom_vline(xintercept = 400, color = 'red') +
  geom_vline(xintercept = 15000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_3@meta.data, aes(x=seurat_obj_3$nFeature_RNA, y=seurat_obj_3$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=30, color='red') + 
  geom_vline(xintercept=300, color='red') +
  geom_vline(xintercept=2700, color='red') 
dev.off()


seurat_obj_3=subset(seurat_obj_3,subset=nFeature_RNA>300 & nFeature_RNA<2700 & MTpercent<30 & nCount_RNA>400 & nCount_RNA<15000)
################################################################################ End Control Healthy MC_2_B #3



################################################################################ Start Control Healthy MC_9 #4
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - Healthy/GSM6250068_MC_9_Colon_661_CD45sort_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - Healthy/GSM6250068_MC_9_Colon_661_CD45sort_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_4 <- CreateSeuratObject(counts = raw_data, project = "Control Healthy MC_9",min.cells=3,min.features=200)
seurat_obj_4[["MTpercent"]]=PercentageFeatureSet(seurat_obj_4,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_4@meta.data, aes(x = seurat_obj_4$nCount_RNA, y = seurat_obj_4$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 30, color = 'red') + 
  geom_vline(xintercept = 280, color = 'red') +
  geom_vline(xintercept = 8000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_4@meta.data, aes(x=seurat_obj_4$nFeature_RNA, y=seurat_obj_4$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=30, color='red') + 
  geom_vline(xintercept=190, color='red') +
  geom_vline(xintercept=2500, color='red') 
dev.off()


seurat_obj_4=subset(seurat_obj_4,subset=nFeature_RNA>190 & nFeature_RNA<2500 & MTpercent<30 & nCount_RNA>280 & nCount_RNA<8000)
################################################################################ End Control Healthy MC_9 #4



################################################################################ Start Control Healthy SIC_13 #5
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - Healthy/GSM6250074_SIC_13_Colon_78_CD45sort_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - Healthy/GSM6250074_SIC_13_Colon_78_CD45sort_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_5 <- CreateSeuratObject(counts = raw_data, project = "Control Healthy SIC_13",min.cells=3,min.features=200)
seurat_obj_5[["MTpercent"]]=PercentageFeatureSet(seurat_obj_5,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_5@meta.data, aes(x = seurat_obj_5$nCount_RNA, y = seurat_obj_5$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 25, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 6000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_5@meta.data, aes(x=seurat_obj_5$nFeature_RNA, y=seurat_obj_5$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=25, color='red') + 
  geom_vline(xintercept=1000, color='red') +
  geom_vline(xintercept=2300, color='red') 
dev.off()


seurat_obj_5=subset(seurat_obj_5,subset=nFeature_RNA>1000 & nFeature_RNA<2300 & MTpercent<25 & nCount_RNA>200 & nCount_RNA<6000)
################################################################################ End Control Healthy SIC_13 #5



################################################################################ Start Control Healthy SIC_186 #6
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - Healthy/GSM6250079_SIC_186_Colon_584_CD45sort_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - Healthy/GSM6250079_SIC_186_Colon_584_CD45sort_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_6 <- CreateSeuratObject(counts = raw_data, project = "Control Healthy SIC_186",min.cells=3,min.features=200)
seurat_obj_6[["MTpercent"]]=PercentageFeatureSet(seurat_obj_6,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_6@meta.data, aes(x = seurat_obj_6$nCount_RNA, y = seurat_obj_6$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 18000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_6@meta.data, aes(x=seurat_obj_6$nFeature_RNA, y=seurat_obj_6$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=400, color='red') +
  geom_vline(xintercept=2300, color='red') 
dev.off()


seurat_obj_6=subset(seurat_obj_6,subset=nFeature_RNA>400 & nFeature_RNA<2300 & MTpercent<20 & nCount_RNA>1000 & nCount_RNA<18000)
################################################################################ End Control Healthy SIC_186 #6



################################################################################ Start Control Healthy SIC_187 #7
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - Healthy/GSM6250080_SIC_187_Colon_584_CD45sort_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - Healthy/GSM6250080_SIC_187_Colon_584_CD45sort_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_7 <- CreateSeuratObject(counts = raw_data, project = "Control Healthy SIC_187",min.cells=3,min.features=200)
seurat_obj_7[["MTpercent"]]=PercentageFeatureSet(seurat_obj_7,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_7@meta.data, aes(x = seurat_obj_7$nCount_RNA, y = seurat_obj_7$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 600, color = 'red') +
  geom_vline(xintercept = 15000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_7@meta.data, aes(x=seurat_obj_7$nFeature_RNA, y=seurat_obj_7$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=300, color='red') +
  geom_vline(xintercept=2000, color='red') 
dev.off()


seurat_obj_7=subset(seurat_obj_7,subset=nFeature_RNA>300 & nFeature_RNA<2000 & MTpercent<20 & nCount_RNA>600 & nCount_RNA<15000)
################################################################################ End Control Healthy SIC_187 #7



################################################################################ Start Control Healthy SIC_188_A #8
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - Healthy/GSM6250081_SIC_188_Colon_585_CD45sort_5p_GEX_A")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - Healthy/GSM6250081_SIC_188_Colon_585_CD45sort_5p_GEX_A/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_8 <- CreateSeuratObject(counts = raw_data, project = "Control Healthy SIC_188_A",min.cells=3,min.features=200)
seurat_obj_8[["MTpercent"]]=PercentageFeatureSet(seurat_obj_8,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_8@meta.data, aes(x = seurat_obj_8$nCount_RNA, y = seurat_obj_8$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 800, color = 'red') +
  geom_vline(xintercept = 15000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_8@meta.data, aes(x=seurat_obj_8$nFeature_RNA, y=seurat_obj_8$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=400, color='red') +
  geom_vline(xintercept=2000, color='red') 
dev.off()


seurat_obj_8=subset(seurat_obj_8,subset=nFeature_RNA>400 & nFeature_RNA<2000 & MTpercent<20 & nCount_RNA>800 & nCount_RNA<15000)
################################################################################ End Control Healthy SIC_188_A #8



################################################################################ Start Control Healthy SIC_188_B #9
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - Healthy/GSM6250082_SIC_188_Colon_585_CD45sort_5p_GEX_B")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - Healthy/GSM6250082_SIC_188_Colon_585_CD45sort_5p_GEX_B/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_9 <- CreateSeuratObject(counts = raw_data, project = "Control Healthy SIC_188_B",min.cells=3,min.features=200)
seurat_obj_9[["MTpercent"]]=PercentageFeatureSet(seurat_obj_9,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_9@meta.data, aes(x = seurat_obj_9$nCount_RNA, y = seurat_obj_9$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 23, color = 'red') + 
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 15000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_9@meta.data, aes(x=seurat_obj_9$nFeature_RNA, y=seurat_obj_9$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=23, color='red') + 
  geom_vline(xintercept=300, color='red') +
  geom_vline(xintercept=2000, color='red') 
dev.off()


seurat_obj_9=subset(seurat_obj_9,subset=nFeature_RNA>300 & nFeature_RNA<2000 & MTpercent<23 & nCount_RNA>700 & nCount_RNA<15000)
################################################################################ End Control Healthy SIC_188_B #9



################################################################################ Start Control Healthy SIC_612_A #10
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - Healthy/GSM6250083_SIC_196_Colon_612_CD45sort_5p_GEX_A")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - Healthy/GSM6250083_SIC_196_Colon_612_CD45sort_5p_GEX_A/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_10 <- CreateSeuratObject(counts = raw_data, project = "Control Healthy SIC_612_A",min.cells=3,min.features=200)
seurat_obj_10[["MTpercent"]]=PercentageFeatureSet(seurat_obj_10,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_10@meta.data, aes(x = seurat_obj_10$nCount_RNA, y = seurat_obj_10$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 25, color = 'red') + 
  geom_vline(xintercept = 500, color = 'red') +
  geom_vline(xintercept = 14000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_10@meta.data, aes(x=seurat_obj_10$nFeature_RNA, y=seurat_obj_10$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=25, color='red') + 
  geom_vline(xintercept=300, color='red') +
  geom_vline(xintercept=1800, color='red') 
dev.off()


seurat_obj_10=subset(seurat_obj_10,subset=nFeature_RNA>300 & nFeature_RNA<1800 & MTpercent<25 & nCount_RNA>500 & nCount_RNA<14000)
################################################################################ End Control Healthy SIC_612_A #10



################################################################################ Start Control Healthy SIC_612_B #11
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - Healthy/GSM6250084_SIC_196_Colon_612_CD45sort_5p_GEX_B")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - Healthy/GSM6250084_SIC_196_Colon_612_CD45sort_5p_GEX_B/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_11 <- CreateSeuratObject(counts = raw_data, project = "Control Healthy SIC_612_B",min.cells=3,min.features=200)
seurat_obj_11[["MTpercent"]]=PercentageFeatureSet(seurat_obj_11,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_11@meta.data, aes(x = seurat_obj_11$nCount_RNA, y = seurat_obj_11$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 25, color = 'red') + 
  geom_vline(xintercept = 400, color = 'red') +
  geom_vline(xintercept = 14000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_11@meta.data, aes(x=seurat_obj_11$nFeature_RNA, y=seurat_obj_11$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=25, color='red') + 
  geom_vline(xintercept=250, color='red') +
  geom_vline(xintercept=1800, color='red') 
dev.off()


seurat_obj_11=subset(seurat_obj_11,subset=nFeature_RNA>150 & nFeature_RNA<1800 & MTpercent<25 & nCount_RNA>400 & nCount_RNA<14000)
################################################################################ End Control Healthy SIC_612_B #11



################################################################################ Start Control - On ICI Therapy SIC_109 #12
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - On ICI Therapy/GSM6250070_SIC_109_Colon_368_CD45pos_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - On ICI Therapy/GSM6250070_SIC_109_Colon_368_CD45pos_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_12 <- CreateSeuratObject(counts = raw_data, project = "Control - On ICI Therapy SIC_109",min.cells=3,min.features=200)
seurat_obj_12[["MTpercent"]]=PercentageFeatureSet(seurat_obj_12,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_12@meta.data, aes(x = seurat_obj_12$nCount_RNA, y = seurat_obj_12$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 47, color = 'red') + 
  geom_vline(xintercept = 500, color = 'red') +
  geom_vline(xintercept = 6000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_12@meta.data, aes(x=seurat_obj_12$nFeature_RNA, y=seurat_obj_12$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=47, color='red') + 
  geom_vline(xintercept=240, color='red') +
  geom_vline(xintercept=1600, color='red') 
dev.off()


seurat_obj_12=subset(seurat_obj_12,subset=nFeature_RNA>240 & nFeature_RNA<1600 & MTpercent<47 & nCount_RNA>500 & nCount_RNA<6000)
################################################################################ End Control - On ICI Therapy SIC_109 #12



################################################################################ Start Control - On ICI Therapy SIC_172 #13
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - On ICI Therapy/GSM6250078_SIC_172_Colon_529_CD45sort_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - On ICI Therapy/GSM6250078_SIC_172_Colon_529_CD45sort_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_13 <- CreateSeuratObject(counts = raw_data, project = "Control - On ICI Therapy SIC_172",min.cells=3,min.features=200)
seurat_obj_13[["MTpercent"]]=PercentageFeatureSet(seurat_obj_13,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_13@meta.data, aes(x = seurat_obj_13$nCount_RNA, y = seurat_obj_13$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 10000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_13@meta.data, aes(x=seurat_obj_13$nFeature_RNA, y=seurat_obj_13$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=600, color='red') +
  geom_vline(xintercept=2900, color='red') 
dev.off()


seurat_obj_13=subset(seurat_obj_13,subset=nFeature_RNA>600 & nFeature_RNA<2900 & MTpercent<20 & nCount_RNA>1000 & nCount_RNA<10000)
################################################################################ End Control - On ICI Therapy SIC_172 #13



################################################################################ Start Control - On ICI Therapy SIC_19 #14
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - On ICI Therapy/GSM6250085_SIC_19_Colon_217_CAP_unselect_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - On ICI Therapy/GSM6250085_SIC_19_Colon_217_CAP_unselect_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_14 <- CreateSeuratObject(counts = raw_data, project = "Control - On ICI Therapy SIC_19",min.cells=3,min.features=200)
seurat_obj_14[["MTpercent"]]=PercentageFeatureSet(seurat_obj_14,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_14@meta.data, aes(x = seurat_obj_14$nCount_RNA, y = seurat_obj_14$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 80, color = 'red') + 
  geom_vline(xintercept = 700, color = 'red') +
  geom_vline(xintercept = 8000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_14@meta.data, aes(x=seurat_obj_14$nFeature_RNA, y=seurat_obj_14$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=80, color='red') + 
  geom_vline(xintercept=240, color='red') +
  geom_vline(xintercept=2600, color='red') 
dev.off()


seurat_obj_14=subset(seurat_obj_14,subset=nFeature_RNA>240 & nFeature_RNA<2600 & MTpercent<80 & nCount_RNA>700 & nCount_RNA<8000)
################################################################################ End Control - On ICI Therapy SIC_19 #14




################################################################################ Start Control - On ICI Therapy SIC_31 #15
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - On ICI Therapy/GSM6250086_SIC_31_Colon_128_CD45sort_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - On ICI Therapy/GSM6250086_SIC_31_Colon_128_CD45sort_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_15 <- CreateSeuratObject(counts = raw_data, project = "Control - On ICI Therapy SIC_31",min.cells=3,min.features=200)
seurat_obj_15[["MTpercent"]]=PercentageFeatureSet(seurat_obj_15,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_15@meta.data, aes(x = seurat_obj_15$nCount_RNA, y = seurat_obj_15$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 18, color = 'red') + 
  geom_vline(xintercept = 900, color = 'red') +
  geom_vline(xintercept = 7000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_15@meta.data, aes(x=seurat_obj_15$nFeature_RNA, y=seurat_obj_15$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=18, color='red') + 
  geom_vline(xintercept=550, color='red') +
  geom_vline(xintercept=2400, color='red') 
dev.off()


seurat_obj_15=subset(seurat_obj_15,subset=nFeature_RNA>550 & nFeature_RNA<2400 & MTpercent<18 & nCount_RNA>900 & nCount_RNA<7000)
################################################################################ End Control - On ICI Therapy SIC_31 #15




################################################################################ Start Control - On ICI Therapy SIC_94 #16
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/Control - On ICI Therapy/GSM6250093_SIC_94_Colon_340_CD45pos_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/Control - On ICI Therapy/GSM6250093_SIC_94_Colon_340_CD45pos_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_16 <- CreateSeuratObject(counts = raw_data, project = "Control - On ICI Therapy SIC_94",min.cells=3,min.features=200)
seurat_obj_16[["MTpercent"]]=PercentageFeatureSet(seurat_obj_16,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_16@meta.data, aes(x = seurat_obj_16$nCount_RNA, y = seurat_obj_16$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 40, color = 'red') + 
  geom_vline(xintercept = 400, color = 'red') +
  geom_vline(xintercept = 3000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_16@meta.data, aes(x=seurat_obj_16$nFeature_RNA, y=seurat_obj_16$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=40, color='red') + 
  geom_vline(xintercept=210, color='red') +
  geom_vline(xintercept=1000, color='red') 
dev.off()


seurat_obj_16=subset(seurat_obj_16,subset=nFeature_RNA>210 & nFeature_RNA<1000 & MTpercent<40 & nCount_RNA>400 & nCount_RNA<3000)
################################################################################ End Control - On ICI Therapy SIC_94 #16



################################################################################ Start irColitis Case SIC_100 #17
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/irColitis Case/GSM6250069_SIC_100_Colon_332_CD45sort_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/irColitis Case/GSM6250069_SIC_100_Colon_332_CD45sort_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_17 <- CreateSeuratObject(counts = raw_data, project = "irColitis Case SIC_100",min.cells=3,min.features=200)
seurat_obj_17[["MTpercent"]]=PercentageFeatureSet(seurat_obj_17,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_17@meta.data, aes(x = seurat_obj_17$nCount_RNA, y = seurat_obj_17$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 500, color = 'red') +
  geom_vline(xintercept = 10000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_17@meta.data, aes(x=seurat_obj_17$nFeature_RNA, y=seurat_obj_17$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=400, color='red') +
  geom_vline(xintercept=3000, color='red') 
dev.off()


seurat_obj_17=subset(seurat_obj_17,subset=nFeature_RNA>400 & nFeature_RNA<3000 & MTpercent<20 & nCount_RNA>500 & nCount_RNA<10000)
################################################################################ End irColitis Case SIC_100 #17




################################################################################ Start irColitis Case SIC_121 #18
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/irColitis Case/GSM6250071_SIC_121_Colon_382_CD45pos_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/irColitis Case/GSM6250071_SIC_121_Colon_382_CD45pos_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_18 <- CreateSeuratObject(counts = raw_data, project = "irColitis Case SIC_121",min.cells=3,min.features=200)
seurat_obj_18[["MTpercent"]]=PercentageFeatureSet(seurat_obj_18,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_18@meta.data, aes(x = seurat_obj_18$nCount_RNA, y = seurat_obj_18$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 25, color = 'red') + 
  geom_vline(xintercept = 300, color = 'red') +
  geom_vline(xintercept = 20000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_18@meta.data, aes(x=seurat_obj_18$nFeature_RNA, y=seurat_obj_18$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=25, color='red') + 
  geom_vline(xintercept=205, color='red') +
  geom_vline(xintercept=4000, color='red') 
dev.off()


seurat_obj_18=subset(seurat_obj_18,subset=nFeature_RNA>205 & nFeature_RNA<4000 & MTpercent<25 & nCount_RNA>300 & nCount_RNA<20000)
################################################################################ End irColitis Case SIC_121 #18



################################################################################ Start irColitis Case SIC_126 #19
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/irColitis Case/GSM6250072_SIC_126_Colon_398_CD45pos_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/irColitis Case/GSM6250072_SIC_126_Colon_398_CD45pos_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_19 <- CreateSeuratObject(counts = raw_data, project = "irColitis Case SIC_126",min.cells=3,min.features=200)
seurat_obj_19[["MTpercent"]]=PercentageFeatureSet(seurat_obj_19,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_19@meta.data, aes(x = seurat_obj_19$nCount_RNA, y = seurat_obj_19$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 350, color = 'red') +
  geom_vline(xintercept = 23000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_19@meta.data, aes(x=seurat_obj_19$nFeature_RNA, y=seurat_obj_19$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=205, color='red') +
  geom_vline(xintercept=3500, color='red') 
dev.off()


seurat_obj_19=subset(seurat_obj_19,subset=nFeature_RNA>205 & nFeature_RNA<3500 & MTpercent<20 & nCount_RNA>350 & nCount_RNA<23000)
################################################################################ End irColitis Case SIC_126 #19




################################################################################ Start irColitis Case SIC_134 #20
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/irColitis Case/GSM6250073_SIC_134_Colon_417_CD45pos_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/irColitis Case/GSM6250073_SIC_134_Colon_417_CD45pos_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_20 <- CreateSeuratObject(counts = raw_data, project = "irColitis Case SIC_134",min.cells=3,min.features=200)
seurat_obj_20[["MTpercent"]]=PercentageFeatureSet(seurat_obj_20,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_20@meta.data, aes(x = seurat_obj_20$nCount_RNA, y = seurat_obj_20$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 18, color = 'red') + 
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 20000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_20@meta.data, aes(x=seurat_obj_20$nFeature_RNA, y=seurat_obj_20$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=18, color='red') + 
  geom_vline(xintercept=205, color='red') +
  geom_vline(xintercept=2500, color='red') 
dev.off()


seurat_obj_20=subset(seurat_obj_20,subset=nFeature_RNA>205 & nFeature_RNA<2500 & MTpercent<18 & nCount_RNA>1000 & nCount_RNA<20000)
################################################################################ End irColitis Case SIC_134 #20




################################################################################ Start irColitis Case SIC_140 #21
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/irColitis Case/GSM6250075_SIC_140_Colon_441_CD45sort_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/irColitis Case/GSM6250075_SIC_140_Colon_441_CD45sort_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_21 <- CreateSeuratObject(counts = raw_data, project = "irColitis Case SIC_140",min.cells=3,min.features=200)
seurat_obj_21[["MTpercent"]]=PercentageFeatureSet(seurat_obj_21,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_21@meta.data, aes(x = seurat_obj_21$nCount_RNA, y = seurat_obj_21$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 18, color = 'red') + 
  geom_vline(xintercept = 1200, color = 'red') +
  geom_vline(xintercept = 23000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_21@meta.data, aes(x=seurat_obj_21$nFeature_RNA, y=seurat_obj_21$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=18, color='red') + 
  geom_vline(xintercept=700, color='red') +
  geom_vline(xintercept=3200, color='red') 
dev.off()


seurat_obj_21=subset(seurat_obj_21,subset=nFeature_RNA>700 & nFeature_RNA<3200 & MTpercent<18 & nCount_RNA>1200 & nCount_RNA<23000)
################################################################################ End irColitis Case SIC_140 #21



################################################################################ Start irColitis Case SIC_141_A #22
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/irColitis Case/GSM6250076_SIC_141_Colon_438_CD45pos_5p_GEX_A")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/irColitis Case/GSM6250076_SIC_141_Colon_438_CD45pos_5p_GEX_A/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_22 <- CreateSeuratObject(counts = raw_data, project = "irColitis Case SIC_141_A",min.cells=3,min.features=200)
seurat_obj_22[["MTpercent"]]=PercentageFeatureSet(seurat_obj_22,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_22@meta.data, aes(x = seurat_obj_22$nCount_RNA, y = seurat_obj_22$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 29, color = 'red') + 
  geom_vline(xintercept = 400, color = 'red') +
  geom_vline(xintercept = 15000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_22@meta.data, aes(x=seurat_obj_22$nFeature_RNA, y=seurat_obj_22$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=29, color='red') + 
  geom_vline(xintercept=200, color='red') +
  geom_vline(xintercept=2800, color='red') 
dev.off()


seurat_obj_22=subset(seurat_obj_22,subset=nFeature_RNA>200 & nFeature_RNA<2800 & MTpercent<29 & nCount_RNA>400 & nCount_RNA<15000)
################################################################################ End irColitis Case SIC_141_A #22



################################################################################ Start irColitis Case SIC_141_B #23
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/irColitis Case/GSM6250077_SIC_141_Colon_438_CD45pos_5p_GEX_B")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/irColitis Case/GSM6250077_SIC_141_Colon_438_CD45pos_5p_GEX_B/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_23 <- CreateSeuratObject(counts = raw_data, project = "irColitis Case SIC_141_B",min.cells=3,min.features=200)
seurat_obj_23[["MTpercent"]]=PercentageFeatureSet(seurat_obj_23,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_23@meta.data, aes(x = seurat_obj_23$nCount_RNA, y = seurat_obj_23$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 29, color = 'red') + 
  geom_vline(xintercept = 400, color = 'red') +
  geom_vline(xintercept = 17000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_23@meta.data, aes(x=seurat_obj_23$nFeature_RNA, y=seurat_obj_23$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=29, color='red') + 
  geom_vline(xintercept=200, color='red') +
  geom_vline(xintercept=2600, color='red') 
dev.off()


seurat_obj_23=subset(seurat_obj_23,subset=nFeature_RNA>200 & nFeature_RNA<2600 & MTpercent<29 & nCount_RNA>400 & nCount_RNA<17000)
################################################################################ End irColitis Case SIC_141_B #23



################################################################################ Start irColitis Case SIC_32_Colon_128 #24
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/irColitis Case/GSM6250087_SIC_32_Colon_128_CAP_unselect_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/irColitis Case/GSM6250087_SIC_32_Colon_128_CAP_unselect_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_24 <- CreateSeuratObject(counts = raw_data, project = "irColitis Case SIC_32_Colon_128",min.cells=3,min.features=200)
seurat_obj_24[["MTpercent"]]=PercentageFeatureSet(seurat_obj_24,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_24@meta.data, aes(x = seurat_obj_24$nCount_RNA, y = seurat_obj_24$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 75, color = 'red') + 
  geom_vline(xintercept = 400, color = 'red') +
  geom_vline(xintercept = 6000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_24@meta.data, aes(x=seurat_obj_24$nFeature_RNA, y=seurat_obj_24$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=75, color='red') + 
  geom_vline(xintercept=200, color='red') +
  geom_vline(xintercept=2000, color='red') 
dev.off()


seurat_obj_24=subset(seurat_obj_24,subset=nFeature_RNA>200 & nFeature_RNA<2000 & MTpercent<75 & nCount_RNA>400 & nCount_RNA<6000)
################################################################################ End irColitis Case SIC_32_Colon_128 #24



################################################################################ Start irColitis Case SIC_32_Colon_178 #25
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/irColitis Case/GSM6250088_SIC_32_Colon_178_CAP_unselect_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/irColitis Case/GSM6250088_SIC_32_Colon_178_CAP_unselect_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_25 <- CreateSeuratObject(counts = raw_data, project = "irColitis Case SIC_32_Colon_178",min.cells=3,min.features=200)
seurat_obj_25[["MTpercent"]]=PercentageFeatureSet(seurat_obj_25,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_25@meta.data, aes(x = seurat_obj_25$nCount_RNA, y = seurat_obj_25$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 75, color = 'red') + 
  geom_vline(xintercept = 400, color = 'red') +
  geom_vline(xintercept = 9000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_25@meta.data, aes(x=seurat_obj_25$nFeature_RNA, y=seurat_obj_25$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=75, color='red') + 
  geom_vline(xintercept=200, color='red') +
  geom_vline(xintercept=2500, color='red') 
dev.off()


seurat_obj_25=subset(seurat_obj_25,subset=nFeature_RNA>200 & nFeature_RNA<2500 & MTpercent<75 & nCount_RNA>400 & nCount_RNA<9000)
################################################################################ End irColitis Case SIC_32_Colon_178 #25




################################################################################ Start irColitis Case SIC_40 #26
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/irColitis Case/GSM6250089_SIC_40_Colon_282_CD45pos_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/irColitis Case/GSM6250089_SIC_40_Colon_282_CD45pos_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_26 <- CreateSeuratObject(counts = raw_data, project = "irColitis Case SIC_40",min.cells=3,min.features=200)
seurat_obj_26[["MTpercent"]]=PercentageFeatureSet(seurat_obj_26,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_26@meta.data, aes(x = seurat_obj_26$nCount_RNA, y = seurat_obj_26$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 30, color = 'red') + 
  geom_vline(xintercept = 300, color = 'red') +
  geom_vline(xintercept = 1800, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_26@meta.data, aes(x=seurat_obj_26$nFeature_RNA, y=seurat_obj_26$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=30, color='red') + 
  geom_vline(xintercept=200, color='red') +
  geom_vline(xintercept=550, color='red') 
dev.off()


seurat_obj_26=subset(seurat_obj_26,subset=nFeature_RNA>200 & nFeature_RNA<550 & MTpercent<30 & nCount_RNA>300 & nCount_RNA<1800)
################################################################################ End irColitis Case SIC_40 #26



################################################################################ Start irColitis Case SIC_71 #27
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/irColitis Case/GSM6250090_SIC_71_Colon_241_CD45sort_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/irColitis Case/GSM6250090_SIC_71_Colon_241_CD45sort_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_27 <- CreateSeuratObject(counts = raw_data, project = "irColitis Case SIC_71",min.cells=3,min.features=200)
seurat_obj_27[["MTpercent"]]=PercentageFeatureSet(seurat_obj_27,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_27@meta.data, aes(x = seurat_obj_27$nCount_RNA, y = seurat_obj_27$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 300, color = 'red') +
  geom_vline(xintercept = 17000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_27@meta.data, aes(x=seurat_obj_27$nFeature_RNA, y=seurat_obj_27$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=400, color='red') +
  geom_vline(xintercept=3000, color='red') 
dev.off()


seurat_obj_27=subset(seurat_obj_27,subset=nFeature_RNA>400 & nFeature_RNA<3000 & MTpercent<20 & nCount_RNA>300 & nCount_RNA<17000)
################################################################################ End irColitis Case SIC_71 #27



################################################################################ Start irColitis Case SIC_76 #28
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/irColitis Case/GSM6250091_SIC_76_Colon_273_CD45pos_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/irColitis Case/GSM6250091_SIC_76_Colon_273_CD45pos_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_28 <- CreateSeuratObject(counts = raw_data, project = "irColitis Case SIC_76",min.cells=3,min.features=200)
seurat_obj_28[["MTpercent"]]=PercentageFeatureSet(seurat_obj_28,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_28@meta.data, aes(x = seurat_obj_28$nCount_RNA, y = seurat_obj_28$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 25, color = 'red') + 
  geom_vline(xintercept = 300, color = 'red') +
  geom_vline(xintercept = 17000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_28@meta.data, aes(x=seurat_obj_28$nFeature_RNA, y=seurat_obj_28$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=25, color='red') + 
  geom_vline(xintercept=200, color='red') +
  geom_vline(xintercept=4000, color='red') 
dev.off()


seurat_obj_28=subset(seurat_obj_28,subset=nFeature_RNA>200 & nFeature_RNA<4000 & MTpercent<25 & nCount_RNA>300 & nCount_RNA<17000)
################################################################################ End irColitis Case SIC_76 #28



################################################################################ Start irColitis Case SIC_89 #29
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/irColitis Case/GSM6250092_SIC_89_Colon_310_CD45pos_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/irColitis Case/GSM6250092_SIC_89_Colon_310_CD45pos_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_29 <- CreateSeuratObject(counts = raw_data, project = "irColitis Case SIC_89",min.cells=3,min.features=200)
seurat_obj_29[["MTpercent"]]=PercentageFeatureSet(seurat_obj_29,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_29@meta.data, aes(x = seurat_obj_29$nCount_RNA, y = seurat_obj_29$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 30, color = 'red') + 
  geom_vline(xintercept = 400, color = 'red') +
  geom_vline(xintercept = 8000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_29@meta.data, aes(x=seurat_obj_29$nFeature_RNA, y=seurat_obj_29$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=30, color='red') + 
  geom_vline(xintercept=200, color='red') +
  geom_vline(xintercept=1700, color='red') 
dev.off()


seurat_obj_29=subset(seurat_obj_29,subset=nFeature_RNA>200 & nFeature_RNA<1700 & MTpercent<30 & nCount_RNA>400 & nCount_RNA<8000)
################################################################################ End irColitis Case SIC_89 #29




################################################################################ Start irColitis Case SIC_97 #30
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/QC/5p/irColitis Case/GSM6250094_SIC_97_Colon_318_CD45pos_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 2/0_Data/GSE206299/GEX/5p/irColitis Case/GSM6250094_SIC_97_Colon_318_CD45pos_5p_GEX/raw_feature_bc_matrix.h5"


raw_data <- Read10X_h5(data_dir)
seurat_obj_30 <- CreateSeuratObject(counts = raw_data, project = "irColitis Case SIC_97",min.cells=3,min.features=200)
seurat_obj_30[["MTpercent"]]=PercentageFeatureSet(seurat_obj_30,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = seurat_obj_30@meta.data, aes(x = seurat_obj_30$nCount_RNA, y = seurat_obj_30$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 30, color = 'red') + 
  geom_vline(xintercept = 350, color = 'red') +
  geom_vline(xintercept = 13000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=seurat_obj_30@meta.data, aes(x=seurat_obj_30$nFeature_RNA, y=seurat_obj_30$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=30, color='red') + 
  geom_vline(xintercept=200, color='red') +
  geom_vline(xintercept=3500, color='red') 
dev.off()


seurat_obj_30=subset(seurat_obj_30,subset=nFeature_RNA>200 & nFeature_RNA<3500 & MTpercent<30 & nCount_RNA>350 & nCount_RNA<13000)
################################################################################ End irColitis Case SIC_97 #30




################################################################################ Start integration
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2")
merged_obj <- merge(seurat_obj_1, y = list(seurat_obj_2, seurat_obj_3, seurat_obj_4, 
                                           seurat_obj_5, seurat_obj_6, seurat_obj_7, 
                                           seurat_obj_8, seurat_obj_9,seurat_obj_10, 
                                           seurat_obj_11, seurat_obj_12,seurat_obj_13,
                                           seurat_obj_14, seurat_obj_15,seurat_obj_16,
                                           seurat_obj_17, seurat_obj_18,seurat_obj_19,
                                           seurat_obj_20, seurat_obj_21,seurat_obj_22,
                                           seurat_obj_23, seurat_obj_24,seurat_obj_25,
                                           seurat_obj_26, seurat_obj_27,seurat_obj_28,
                                           seurat_obj_29, seurat_obj_30))


merged_obj=NormalizeData(merged_obj,normalization.method = "LogNormalize",scale.factor = 10000)


merged_obj=FindVariableFeatures(merged_obj,selection.method = "vst",nfeatures = 2000)


merged_obj = ScaleData(merged_obj,features = rownames(merged_obj))


merged_obj = RunPCA(merged_obj)

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


# 3
# saveRDS(file = "merged_obj1",merged_obj1)
# The Seurat object obtained after FindNeighbors and FindClusters


merged_obj1=RunUMAP(merged_obj1,dims = 1:30, reduction = "integrated.cca")
png(filename = "UMAP1.png",width = 9000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1,label = TRUE)
dev.off()


png(filename = "FeaturePlot.png", width = 9000, height = 4000, units = "px", res = 600)
FeaturePlot(merged_obj1, features = c("CD3D", "CD3E", "CD3G", "LYZ", "CD79A", "CD19"))
dev.off()


# Remove clusters 0 2 5 6 8 9 10
merged_obj1 <- subset(merged_obj1, subset = seurat_clusters %in% c(0, 2, 5, 6, 8, 9, 10), invert = TRUE)
png(filename = "UMAP2.png",width = 7000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1,label = TRUE)
dev.off()


# Remove cells with an x-axis value greater than 1.
umap_coord <- merged_obj1@reductions[["umap"]]@cell.embeddings
cells_to_keep <- rownames(umap_coord[umap_coord[, "umap_1"] <= 1, ])
merged_obj1 <- subset(merged_obj1, cells = cells_to_keep)
png(filename = "UMAP3.png", width = 7000, height = 4000, units = "px", res = 600)
DimPlot(merged_obj1, label = TRUE)
dev.off()


# This code changes the label of cluster 7 to cluster 2
#  and then generates and saves a UMAP plot of the modified dataset.
cluster_ids <- as.character(Idents(merged_obj1))
cluster_ids[cluster_ids == "7"] <- "2" 

Idents(merged_obj1) <- cluster_ids
png(filename = "UMAP4.png", width = 7000, height = 4000, units = "px", res = 600)
DimPlot(merged_obj1, label = TRUE)
dev.off()


# This code reorders the cluster identities in the Seurat object (merged_obj1) according to the specified order ("1", "2", "3", "4")
desired_order <- c("1", "2", "3", "4")
Idents(merged_obj1) <- factor(Idents(merged_obj1), levels = desired_order)
png(filename = "UMAP5.png", width = 7000, height = 4000, units = "px", res = 600)
DimPlot(merged_obj1, label = TRUE)
dev.off()


png(filename = "FeaturePlot2.png", width = 5000, height = 4000, units = "px", res = 600)
FeaturePlot(merged_obj1, features = c("CD3D", "CD3E", "CD3G", "LYZ", "CD79A", "CD19"))
dev.off()


# 4
# saveRDS(file = "merged_obj1",merged_obj1)
# The Seurat object obtained after UMAP and FeaturePlot

################################################################################ End UMAP

merged_obj2 = merged_obj1

################################################################################ Start Extracting and saving Seurat objects for each sample

samples <- unique(merged_obj2$orig.ident)

for (i in seq_along(samples)) {
  sample_name <- samples[i]
  seurat_obj <- subset(merged_obj2, subset = orig.ident == sample_name)
  assign(paste0("seurat_obj_", i), seurat_obj)
}


# 5
# GSE206299_seurat_objs
seurat_objs <- list(seurat_obj_1, seurat_obj_2, seurat_obj_3, seurat_obj_4, seurat_obj_5, 
                    seurat_obj_6, seurat_obj_7, seurat_obj_8, seurat_obj_9, seurat_obj_10, 
                    seurat_obj_11, seurat_obj_12, seurat_obj_13, seurat_obj_14, seurat_obj_15, 
                    seurat_obj_16, seurat_obj_17, seurat_obj_18, seurat_obj_19, seurat_obj_20, 
                    seurat_obj_21, seurat_obj_22, seurat_obj_23, seurat_obj_24, seurat_obj_25, 
                    seurat_obj_26, seurat_obj_27, seurat_obj_28, seurat_obj_29, seurat_obj_30)

for (i in 1:length(seurat_objs)) {
  saveRDS(seurat_objs[[i]], file = paste0("seurat_obj_", i, ".rds"))
}

################################################################################ End Extracting and saving Seurat objects for each sample