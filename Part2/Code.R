################################################################################ Start library
library(Seurat)
library(SeuratObject)
library(ggplot2)

install.packages("BiocManager")
BiocManager::install("rhdf5")
library(rhdf5)
################################################################################ End library



################################################################################ Start Control Healthy MC_1 #1
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250065_MC_1_Colon_555_CD45pos_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250065_MC_1_Colon_555_CD45pos_5p_GEX/raw_feature_bc_matrix.h5"


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250066_MC_2_Colon_576_CD45sort_5p_GEX_A")
data_dir <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250066_MC_2_Colon_576_CD45sort_5p_GEX_A/raw_feature_bc_matrix.h5"


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250067_MC_2_Colon_576_CD45sort_5p_GEX_B")
data_dir <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250067_MC_2_Colon_576_CD45sort_5p_GEX_B/raw_feature_bc_matrix.h5"


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250068_MC_9_Colon_661_CD45sort_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250068_MC_9_Colon_661_CD45sort_5p_GEX/raw_feature_bc_matrix.h5"


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250074_SIC_13_Colon_78_CD45sort_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250074_SIC_13_Colon_78_CD45sort_5p_GEX/raw_feature_bc_matrix.h5"


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250079_SIC_186_Colon_584_CD45sort_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250079_SIC_186_Colon_584_CD45sort_5p_GEX/raw_feature_bc_matrix.h5"


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250080_SIC_187_Colon_584_CD45sort_5p_GEX")
data_dir <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250080_SIC_187_Colon_584_CD45sort_5p_GEX/raw_feature_bc_matrix.h5"


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250081_SIC_188_Colon_585_CD45sort_5p_GEX_A")
data_dir <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250081_SIC_188_Colon_585_CD45sort_5p_GEX_A/raw_feature_bc_matrix.h5"


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250082_SIC_188_Colon_585_CD45sort_5p_GEX_B")
data_dir <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250082_SIC_188_Colon_585_CD45sort_5p_GEX_B/raw_feature_bc_matrix.h5"


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250083_SIC_196_Colon_612_CD45sort_5p_GEX_A")
data_dir <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250083_SIC_196_Colon_612_CD45sort_5p_GEX_A/raw_feature_bc_matrix.h5"


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250084_SIC_196_Colon_612_CD45sort_5p_GEX_B")
data_dir <- "C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part2/Data/5p/Control - Healthy/GSM6250084_SIC_196_Colon_612_CD45sort_5p_GEX_B/raw_feature_bc_matrix.h5"


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