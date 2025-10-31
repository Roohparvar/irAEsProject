################################################################################ Start libraries
library(Seurat)
library(SeuratObject)
library(ggplot2)
################################################################################ End libraries



################################################################################ Start Normal Control CT1 #1
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/Normal Control/GSM4288840_CT1-CD3-genes-barcodes-matrix")


srobj_1=CreateSeuratObject(CountMatrix,project ="Normal Control CT1",min.cells=3,min.features=200)
srobj_1[["MTpercent"]]=PercentageFeatureSet(srobj_1,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/Normal Control/GSM4288840_CT1-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_1@meta.data, aes(x = srobj_1$nCount_RNA, y = srobj_1$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 25, color = 'red') + 
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 25000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_1@meta.data, aes(x=srobj_1$nFeature_RNA, y=srobj_1$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=25, color='red') + 
  geom_vline(xintercept=600, color='red') +
  geom_vline(xintercept=5000, color='red') 
dev.off()


srobj_1=subset(srobj_1,subset=nFeature_RNA>600 & nFeature_RNA<5000 & MTpercent<20 & nCount_RNA>1000 & nCount_RNA<25000)
################################################################################ End Normal Control CT1 #1



################################################################################ Start Normal Control CT2 #2
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/Normal Control/GSM4288841_CT2-CD3-genes-barcodes-matrix")


srobj_2=CreateSeuratObject(CountMatrix,project ="Normal Control CT2",min.cells=3,min.features=200)
srobj_2[["MTpercent"]]=PercentageFeatureSet(srobj_2,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/Normal Control/GSM4288841_CT2-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_2@meta.data, aes(x = srobj_2$nCount_RNA, y = srobj_2$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1500, color = 'red') +
  geom_vline(xintercept = 21000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_2@meta.data, aes(x=srobj_2$nFeature_RNA, y=srobj_2$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=600, color='red') +
  geom_vline(xintercept=5000, color='red') 
dev.off()


srobj_2=subset(srobj_2,subset=nFeature_RNA>600 & nFeature_RNA<5000 & MTpercent<20 & nCount_RNA>1500 & nCount_RNA<21000)
################################################################################ End Normal Control CT2 #2



################################################################################ Start Normal Control CT3 #3
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/Normal Control/GSM4288842_CT3-CD3-genes-barcodes-matrix")


srobj_3=CreateSeuratObject(CountMatrix,project ="Normal Control CT3",min.cells=3,min.features=200)
srobj_3[["MTpercent"]]=PercentageFeatureSet(srobj_3,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/Normal Control/GSM4288842_CT3-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_3@meta.data, aes(x = srobj_3$nCount_RNA, y = srobj_3$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1800, color = 'red') +
  geom_vline(xintercept = 20000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_3@meta.data, aes(x=srobj_3$nFeature_RNA, y=srobj_3$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=5000, color='red') 
dev.off()


srobj_3=subset(srobj_3,subset=nFeature_RNA>800 & nFeature_RNA<5000 & MTpercent<20 & nCount_RNA>1800 & nCount_RNA<20000)
################################################################################ End Normal Control CT3 #3



################################################################################ Start Normal Control CT4 #4
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/Normal Control/GSM4288843_CT4-CD3-genes-barcodes-matrix")


srobj_4=CreateSeuratObject(CountMatrix,project ="Normal Control CT4",min.cells=3,min.features=200)
srobj_4[["MTpercent"]]=PercentageFeatureSet(srobj_4,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/Normal Control/GSM4288843_CT4-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_4@meta.data, aes(x = srobj_4$nCount_RNA, y = srobj_4$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 18, color = 'red') + 
  geom_vline(xintercept = 1800, color = 'red') +
  geom_vline(xintercept = 18000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_4@meta.data, aes(x=srobj_4$nFeature_RNA, y=srobj_4$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=18, color='red') + 
  geom_vline(xintercept=900, color='red') +
  geom_vline(xintercept=4500, color='red') 
dev.off()


srobj_4=subset(srobj_4,subset=nFeature_RNA>900 & nFeature_RNA<4500 & MTpercent<18 & nCount_RNA>1800 & nCount_RNA<18000)
################################################################################ End Normal Control CT4 #4



################################################################################ Start Normal Control CT5 #5
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/Normal Control/GSM4288844_CT5-CD3-genes-barcodes-matrix")


srobj_5=CreateSeuratObject(CountMatrix,project ="Normal Control CT5",min.cells=3,min.features=200)
srobj_5[["MTpercent"]]=PercentageFeatureSet(srobj_5,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/Normal Control/GSM4288844_CT5-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_5@meta.data, aes(x = srobj_5$nCount_RNA, y = srobj_5$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 25, color = 'red') + 
  geom_vline(xintercept = 1800, color = 'red') +
  geom_vline(xintercept = 25000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_5@meta.data, aes(x=srobj_5$nFeature_RNA, y=srobj_5$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=25, color='red') + 
  geom_vline(xintercept=900, color='red') +
  geom_vline(xintercept=4800, color='red') 
dev.off()


srobj_5=subset(srobj_5,subset=nFeature_RNA>900 & nFeature_RNA<4800 & MTpercent<25 & nCount_RNA>1800 & nCount_RNA<25000)
################################################################################ End Normal Control CT5 #5



################################################################################ Start Normal Control CT6 #6
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/Normal Control/GSM4288845_CT6-CD3-genes-barcodes-matrix")


srobj_6=CreateSeuratObject(CountMatrix,project ="Normal Control CT6",min.cells=3,min.features=200)
srobj_6[["MTpercent"]]=PercentageFeatureSet(srobj_6,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/Normal Control/GSM4288845_CT6-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_6@meta.data, aes(x = srobj_6$nCount_RNA, y = srobj_6$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 23, color = 'red') + 
  geom_vline(xintercept = 1500, color = 'red') +
  geom_vline(xintercept = 18000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_6@meta.data, aes(x=srobj_6$nFeature_RNA, y=srobj_6$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=23, color='red') + 
  geom_vline(xintercept=700, color='red') +
  geom_vline(xintercept=4000, color='red') 
dev.off()


srobj_6=subset(srobj_6,subset=nFeature_RNA>700 & nFeature_RNA<4000 & MTpercent<23 & nCount_RNA>1500 & nCount_RNA<18000)
################################################################################ End Normal Control CT6 #6



################################################################################ Start Normal Control CT7 #7
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/Normal Control/GSM4288846_CT7-CD3-genes-barcodes-matrix")


srobj_7=CreateSeuratObject(CountMatrix,project ="Normal Control CT7",min.cells=3,min.features=200)
srobj_7[["MTpercent"]]=PercentageFeatureSet(srobj_7,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/Normal Control/GSM4288846_CT7-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_7@meta.data, aes(x = srobj_7$nCount_RNA, y = srobj_7$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 21, color = 'red') + 
  geom_vline(xintercept = 1400, color = 'red') +
  geom_vline(xintercept = 19000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_7@meta.data, aes(x=srobj_7$nFeature_RNA, y=srobj_7$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=21, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=4000, color='red') 
dev.off()


srobj_7=subset(srobj_7,subset=nFeature_RNA>800 & nFeature_RNA<4000 & MTpercent<21 & nCount_RNA>1400 & nCount_RNA<19000)
################################################################################ End Normal Control CT7 #7



################################################################################ Start Normal Control CT8 #8
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/Normal Control/GSM4288847_CT8-CD3-genes-barcodes-matrix")


srobj_8=CreateSeuratObject(CountMatrix,project ="Normal Control CT8",min.cells=3,min.features=200)
srobj_8[["MTpercent"]]=PercentageFeatureSet(srobj_8,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/Normal Control/GSM4288847_CT8-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_8@meta.data, aes(x = srobj_8$nCount_RNA, y = srobj_8$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1400, color = 'red') +
  geom_vline(xintercept = 17000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_8@meta.data, aes(x=srobj_8$nFeature_RNA, y=srobj_8$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=4000, color='red') 
dev.off()


srobj_8=subset(srobj_8,subset=nFeature_RNA>800 & nFeature_RNA<4000 & MTpercent<20 & nCount_RNA>1400 & nCount_RNA<17000)
################################################################################ End Normal Control CT8 #8



################################################################################ Start +CPI no colitis NC1 #9
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/+CPI no colitis/GSM4288834_NC1-CD3-genes-barcodes-matrix")


srobj_9=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC1",min.cells=3,min.features=200)
srobj_9[["MTpercent"]]=PercentageFeatureSet(srobj_9,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/+CPI no colitis/GSM4288834_NC1-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_9@meta.data, aes(x = srobj_9$nCount_RNA, y = srobj_9$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 15, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 21000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_9@meta.data, aes(x=srobj_9$nFeature_RNA, y=srobj_9$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=15, color='red') + 
  geom_vline(xintercept=900, color='red') +
  geom_vline(xintercept=6000, color='red') 
dev.off()


srobj_9=subset(srobj_9,subset=nFeature_RNA>900 & nFeature_RNA<6000 & MTpercent<15 & nCount_RNA>2000 & nCount_RNA<21000)
################################################################################ End +CPI no colitis NC1 #9



################################################################################ Start +CPI no colitis NC2 #10
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/+CPI no colitis/GSM4288835_NC2-CD3-genes-barcodes-matrix")


srobj_10=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC2",min.cells=3,min.features=200)
srobj_10[["MTpercent"]]=PercentageFeatureSet(srobj_10,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/+CPI no colitis/GSM4288835_NC2-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_10@meta.data, aes(x = srobj_10$nCount_RNA, y = srobj_10$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 19, color = 'red') + 
  geom_vline(xintercept = 1600, color = 'red') +
  geom_vline(xintercept = 18000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_10@meta.data, aes(x=srobj_10$nFeature_RNA, y=srobj_10$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=19, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=4500, color='red') 
dev.off()


srobj_10=subset(srobj_10,subset=nFeature_RNA>800 & nFeature_RNA<4500 & MTpercent<19 & nCount_RNA>1600 & nCount_RNA<18000)
################################################################################ End +CPI no colitis NC2 #10



################################################################################ Start +CPI no colitis NC3 #11
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/+CPI no colitis/GSM4288836_NC3-CD3-genes-barcodes-matrix")


srobj_11=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC3",min.cells=3,min.features=200)
srobj_11[["MTpercent"]]=PercentageFeatureSet(srobj_11,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/+CPI no colitis/GSM4288836_NC3-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_11@meta.data, aes(x = srobj_11$nCount_RNA, y = srobj_11$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1600, color = 'red') +
  geom_vline(xintercept = 18000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_11@meta.data, aes(x=srobj_11$nFeature_RNA, y=srobj_11$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=4200, color='red') 
dev.off()


srobj_11=subset(srobj_11,subset=nFeature_RNA>800 & nFeature_RNA<4200 & MTpercent<20 & nCount_RNA>1600 & nCount_RNA<18000)
################################################################################ End +CPI no colitis NC3 #11



################################################################################ Start +CPI no colitis NC4 #12
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/+CPI no colitis/GSM4288837_NC4-CD3-genes-barcodes-matrix")


srobj_12=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC4",min.cells=3,min.features=200)
srobj_12[["MTpercent"]]=PercentageFeatureSet(srobj_12,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/+CPI no colitis/GSM4288837_NC4-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_12@meta.data, aes(x = srobj_12$nCount_RNA, y = srobj_12$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1600, color = 'red') +
  geom_vline(xintercept = 18000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_12@meta.data, aes(x=srobj_12$nFeature_RNA, y=srobj_12$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=4200, color='red') 
dev.off()


srobj_12=subset(srobj_12,subset=nFeature_RNA>800 & nFeature_RNA<4200 & MTpercent<20 & nCount_RNA>1600 & nCount_RNA<18000)
################################################################################ End +CPI no colitis NC4 #12



################################################################################ Start +CPI no colitis NC5 #13
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/+CPI no colitis/GSM4288838_NC5-CD3-genes-barcodes-matrix")


srobj_13=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC5",min.cells=3,min.features=200)
srobj_13[["MTpercent"]]=PercentageFeatureSet(srobj_13,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/+CPI no colitis/GSM4288838_NC5-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_13@meta.data, aes(x = srobj_13$nCount_RNA, y = srobj_13$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 23, color = 'red') + 
  geom_vline(xintercept = 1600, color = 'red') +
  geom_vline(xintercept = 22000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_13@meta.data, aes(x=srobj_13$nFeature_RNA, y=srobj_13$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=23, color='red') + 
  geom_vline(xintercept=900, color='red') +
  geom_vline(xintercept=5000, color='red') 
dev.off()


srobj_13=subset(srobj_13,subset=nFeature_RNA>900 & nFeature_RNA<5000 & MTpercent<23 & nCount_RNA>1600 & nCount_RNA<22000)
################################################################################ End +CPI no colitis NC5 #13




################################################################################ Start +CPI no colitis NC6 #14
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/+CPI no colitis/GSM4288839_NC6-CD3-genes-barcodes-matrix")


srobj_14=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC6",min.cells=3,min.features=200)
srobj_14[["MTpercent"]]=PercentageFeatureSet(srobj_14,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/+CPI no colitis/GSM4288839_NC6-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_14@meta.data, aes(x = srobj_14$nCount_RNA, y = srobj_14$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 21, color = 'red') + 
  geom_vline(xintercept = 1400, color = 'red') +
  geom_vline(xintercept = 21000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_14@meta.data, aes(x=srobj_14$nFeature_RNA, y=srobj_14$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=21, color='red') + 
  geom_vline(xintercept=600, color='red') +
  geom_vline(xintercept=5000, color='red') 
dev.off()


srobj_14=subset(srobj_14,subset=nFeature_RNA>600 & nFeature_RNA<5000 & MTpercent<21 & nCount_RNA>1400 & nCount_RNA<21000)
################################################################################ End +CPI no colitis NC6 #14




################################################################################ Start +CPI colitis C1 #15
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/+CPI colitis/GSM4288826_C1-CD3-genes-barcodes-matrix")


srobj_15=CreateSeuratObject(CountMatrix,project ="+CPI colitis C1",min.cells=3,min.features=200)
srobj_15[["MTpercent"]]=PercentageFeatureSet(srobj_15,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/+CPI colitis/GSM4288826_C1-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_15@meta.data, aes(x = srobj_15$nCount_RNA, y = srobj_15$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 25, color = 'red') + 
  geom_vline(xintercept = 1600, color = 'red') +
  geom_vline(xintercept = 30000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_15@meta.data, aes(x=srobj_15$nFeature_RNA, y=srobj_15$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=25, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=6200, color='red') 
dev.off()


srobj_15=subset(srobj_15,subset=nFeature_RNA>800 & nFeature_RNA<6200 & MTpercent<25 & nCount_RNA>1600 & nCount_RNA<30000)
################################################################################ End +CPI colitis C1 #15



################################################################################ Start +CPI colitis C2 #16
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/+CPI colitis/GSM4288827_C2-CD3-genes-barcodes-matrix")


srobj_16=CreateSeuratObject(CountMatrix,project ="+CPI colitis C2",min.cells=3,min.features=200)
srobj_16[["MTpercent"]]=PercentageFeatureSet(srobj_16,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/+CPI colitis/GSM4288827_C2-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_16@meta.data, aes(x = srobj_16$nCount_RNA, y = srobj_16$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1600, color = 'red') +
  geom_vline(xintercept = 37000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_16@meta.data, aes(x=srobj_16$nFeature_RNA, y=srobj_16$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=900, color='red') +
  geom_vline(xintercept=6200, color='red') 
dev.off()


srobj_16=subset(srobj_16,subset=nFeature_RNA>900 & nFeature_RNA<6200 & MTpercent<20 & nCount_RNA>1600 & nCount_RNA<37000)
################################################################################ End +CPI colitis C2 #16



################################################################################ Start +CPI colitis C3 #17
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/+CPI colitis/GSM4288828_C3-CD3-genes-barcodes-matrix")


srobj_17=CreateSeuratObject(CountMatrix,project ="+CPI colitis C3",min.cells=3,min.features=200)
srobj_17[["MTpercent"]]=PercentageFeatureSet(srobj_17,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/+CPI colitis/GSM4288828_C3-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_17@meta.data, aes(x = srobj_17$nCount_RNA, y = srobj_17$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1400, color = 'red') +
  geom_vline(xintercept = 37000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_17@meta.data, aes(x=srobj_17$nFeature_RNA, y=srobj_17$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=6200, color='red') 
dev.off()


srobj_17=subset(srobj_17,subset=nFeature_RNA>800 & nFeature_RNA<6200 & MTpercent<20 & nCount_RNA>1400 & nCount_RNA<37000)
################################################################################ End +CPI colitis C3 #17



################################################################################ Start +CPI colitis C4 #18
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/+CPI colitis/GSM4288829_C4-CD3-genes-barcodes-matrix")


srobj_18=CreateSeuratObject(CountMatrix,project ="+CPI colitis C4",min.cells=3,min.features=200)
srobj_18[["MTpercent"]]=PercentageFeatureSet(srobj_18,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/+CPI colitis/GSM4288829_C4-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_18@meta.data, aes(x = srobj_18$nCount_RNA, y = srobj_18$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1400, color = 'red') +
  geom_vline(xintercept = 37000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_18@meta.data, aes(x=srobj_18$nFeature_RNA, y=srobj_18$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=6200, color='red') 
dev.off()


srobj_18=subset(srobj_18,subset=nFeature_RNA>800 & nFeature_RNA<6200 & MTpercent<20 & nCount_RNA>1400 & nCount_RNA<37000)
################################################################################ End +CPI colitis C4 #18



################################################################################ Start +CPI colitis C5 #19
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/+CPI colitis/GSM4288830_C5-CD3-genes-barcodes-matrix")


srobj_19=CreateSeuratObject(CountMatrix,project ="+CPI colitis C5",min.cells=3,min.features=200)
srobj_19[["MTpercent"]]=PercentageFeatureSet(srobj_19,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/+CPI colitis/GSM4288830_C5-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_19@meta.data, aes(x = srobj_19$nCount_RNA, y = srobj_19$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 17, color = 'red') + 
  geom_vline(xintercept = 1400, color = 'red') +
  geom_vline(xintercept = 37000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_19@meta.data, aes(x=srobj_19$nFeature_RNA, y=srobj_19$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=17, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=5600, color='red') 
dev.off()


srobj_19=subset(srobj_19,subset=nFeature_RNA>800 & nFeature_RNA<5600 & MTpercent<17 & nCount_RNA>1400 & nCount_RNA<37000)
################################################################################ End +CPI colitis C5 #19



################################################################################ Start +CPI colitis C6 #20
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/+CPI colitis/GSM4288831_C6-CD3-genes-barcodes-matrix")


srobj_20=CreateSeuratObject(CountMatrix,project ="+CPI colitis C6",min.cells=3,min.features=200)
srobj_20[["MTpercent"]]=PercentageFeatureSet(srobj_20,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/+CPI colitis/GSM4288831_C6-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_20@meta.data, aes(x = srobj_20$nCount_RNA, y = srobj_20$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1400, color = 'red') +
  geom_vline(xintercept = 37000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_20@meta.data, aes(x=srobj_20$nFeature_RNA, y=srobj_20$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=5600, color='red') 
dev.off()


srobj_20=subset(srobj_20,subset=nFeature_RNA>800 & nFeature_RNA<5600 & MTpercent<20 & nCount_RNA>1400 & nCount_RNA<37000)
################################################################################ End +CPI colitis C6 #20



################################################################################ Start +CPI colitis C7 #21
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/+CPI colitis/GSM4288832_C7-CD3-genes-barcodes-matrix")


srobj_21=CreateSeuratObject(CountMatrix,project ="+CPI colitis C7",min.cells=3,min.features=200)
srobj_21[["MTpercent"]]=PercentageFeatureSet(srobj_21,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/+CPI colitis/GSM4288832_C7-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_21@meta.data, aes(x = srobj_21$nCount_RNA, y = srobj_21$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 18, color = 'red') + 
  geom_vline(xintercept = 1600, color = 'red') +
  geom_vline(xintercept = 43000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_21@meta.data, aes(x=srobj_21$nFeature_RNA, y=srobj_21$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=18, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=6500, color='red') 
dev.off()


srobj_21=subset(srobj_21,subset=nFeature_RNA>800 & nFeature_RNA<6500 & MTpercent<18 & nCount_RNA>1600 & nCount_RNA<43000)
################################################################################ End +CPI colitis C7 #21



################################################################################ Start +CPI colitis C8 #22
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part 1/0_Data/GEX_Data/CD3/+CPI colitis/GSM4288833_C8-CD3-genes-barcodes-matrix")


srobj_22=CreateSeuratObject(CountMatrix,project ="+CPI colitis C8",min.cells=3,min.features=200)
srobj_22[["MTpercent"]]=PercentageFeatureSet(srobj_22,pattern = "^MT-")


setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1/QC/CD3/+CPI colitis/GSM4288833_C8-CD3-genes-barcodes-matrix")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_22@meta.data, aes(x = srobj_22$nCount_RNA, y = srobj_22$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 25, color = 'red') + 
  geom_vline(xintercept = 1100, color = 'red') +
  geom_vline(xintercept = 47000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_22@meta.data, aes(x=srobj_22$nFeature_RNA, y=srobj_22$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=25, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=6500, color='red') 
dev.off()


srobj_22=subset(srobj_22,subset=nFeature_RNA>800 & nFeature_RNA<6500 & MTpercent<25 & nCount_RNA>1100 & nCount_RNA<47000)
################################################################################ End +CPI colitis C8 #22



################################################################################ Start integration and UMAP
setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part 1")
merged_obj <- merge(srobj_1, y = list(srobj_2, srobj_3, srobj_4, 
                                      srobj_5, srobj_6, srobj_7, 
                                      srobj_8, srobj_9,srobj_10, 
                                      srobj_11, srobj_12,srobj_13,
                                      srobj_14, srobj_15,srobj_16,
                                      srobj_17, srobj_18,srobj_19,
                                      srobj_20, srobj_21,srobj_22))


merged_obj=NormalizeData(merged_obj,normalization.method = "LogNormalize",scale.factor = 10000)


merged_obj=FindVariableFeatures(merged_obj,selection.method = "vst",nfeatures = 2000)


merged_obj = ScaleData(merged_obj,features = rownames(merged_obj))


merged_obj = RunPCA(merged_obj)

# 1
# saveRDS(file = "merged_obj",merged_obj)
# The Seurat object obtained after RunPCA and before IntegrateLayers
 
total_barcodes <- length(colnames(merged_obj))
unique_barcodes <- length(unique(colnames(merged_obj)))
 
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


merged_obj1=RunUMAP(merged_obj1,dims = 1:30, reduction = "integrated.cca")
png(filename = "1_First UMAP.png", width = 4000,height = 3000, units ="px", res = 600 )
DimPlot(merged_obj1,label = TRUE)
dev.off()


png(filename = "2_Feature Plot.png", width = 5000, height = 6000, units = "px", res = 600)
FeaturePlot(merged_obj1, features = c("CD3D", "CD3E", "CD3G", "LYZ", "CD79A", "CD19"))
dev.off()
################################################################################ End UMAP

merged_obj2 = merged_obj1


merged_obj2$orig.ident <- dplyr::case_when(
  grepl("Normal Control CT", merged_obj2$orig.ident) ~ sub("Normal Control ", "Healthy ", merged_obj2$orig.ident),
  grepl("\\+CPI no colitis NC", merged_obj2$orig.ident) ~ sub("\\+CPI no colitis", "CPI_Control", merged_obj2$orig.ident),
  grepl("\\+CPI colitis C", merged_obj2$orig.ident) ~ sub("\\+CPI colitis", "CPI_Colitis", merged_obj2$orig.ident),
  TRUE ~ merged_obj2$orig.ident
)
merged_obj2 <- RenameCells(merged_obj2, new.names = paste0("Adrienne/", merged_obj2$orig.ident, "/", colnames(merged_obj2)))



cell_counts <- as.data.frame(table(merged_obj2$orig.ident))
colnames(cell_counts) <- c("Sample", "Number_of_Cells")
print(cell_counts, row.names = FALSE)

# CPI_Colitis C1 : 3394
# CPI_Colitis C2 : 3889
# CPI_Colitis C3: 4128
# CPI_Colitis C4 : 3445
# CPI_Colitis C5 : 3892
# CPI_Colitis C6 : 3011
# CPI_Colitis C7 : 3465
# CPI_Colitis C8 : 3517
# CPI_Control NC1 : 2918
# CPI_Control NC2 : 3661
# CPI_Control NC3 : 2568
# CPI_Control NC4 : 3465
# CPI_Control NC5 : 3222
# CPI_Control NC6 : 3995
# Healthy CT1 : 3708
# Healthy CT2 : 2637
# Healthy CT3 : 3652
# Healthy CT4 : 2479
# Healthy CT5 : 2690
# Healthy CT6 : 4598
# Healthy CT7 : 3408
# Healthy CT8 : 3148

################################################################################ Start Extracting and saving Seurat objects for each sample
setwd("C:/Esmaeil/irAEsProject/Backup/Part 1/3_The Seurat objects per sample")
samples <- unique(merged_obj2$orig.ident)

for (i in seq_along(samples)) {
  sample_name <- samples[i]
  seurat_obj <- subset(merged_obj2, subset = orig.ident == sample_name)
  assign(paste0("srobj_", i), seurat_obj)
}


# 3
# 3_The Seurat objects per sample
seurat_objs <- list(srobj_1, srobj_2, srobj_3, srobj_4, srobj_5, 
                    srobj_6, srobj_7, srobj_8, srobj_9, srobj_10, 
                    srobj_11, srobj_12, srobj_13, srobj_14, srobj_15, 
                    srobj_16, srobj_17, srobj_18, srobj_19, srobj_20, 
                    srobj_21, srobj_22)

# Loop to save each Seurat object
for (i in 1:length(seurat_objs)) {
  saveRDS(seurat_objs[[i]], file = paste0("srobj_", i, ".rds"))
}

################################################################################ End Extracting and saving Seurat objects for each sample