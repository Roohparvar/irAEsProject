################################################################################ Start library
library(Seurat)
library(SeuratObject)
library(ggplot2)
################################################################################ End library



################################################################################ Start Normal Control CT1 #1
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/Normal Control/GSM4288840_CT1-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/Normal Control/GSM4288840_CT1-CD3-genes-barcodes-matrix")


srobj_1=CreateSeuratObject(CountMatrix,project ="Normal Control CT1",min.cells=3,min.features=200)
srobj_1[["MTpercent"]]=PercentageFeatureSet(srobj_1,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/GSE144469/Data/CD3/Normal Control/GSM4288841_CT2-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/GSE144469/Data/CD3/Normal Control/GSM4288841_CT2-CD3-genes-barcodes-matrix")


srobj_2=CreateSeuratObject(CountMatrix,project ="Normal Control CT2",min.cells=3,min.features=200)
srobj_2[["MTpercent"]]=PercentageFeatureSet(srobj_2,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/Normal Control/GSM4288842_CT3-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/Normal Control/GSM4288842_CT3-CD3-genes-barcodes-matrix")


srobj_3=CreateSeuratObject(CountMatrix,project ="Normal Control CT3",min.cells=3,min.features=200)
srobj_3[["MTpercent"]]=PercentageFeatureSet(srobj_3,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/Normal Control/GSM4288843_CT4-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/Normal Control/GSM4288843_CT4-CD3-genes-barcodes-matrix")


srobj_4=CreateSeuratObject(CountMatrix,project ="Normal Control CT4",min.cells=3,min.features=200)
srobj_4[["MTpercent"]]=PercentageFeatureSet(srobj_4,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/Normal Control/GSM4288844_CT5-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/Normal Control/GSM4288844_CT5-CD3-genes-barcodes-matrix")


srobj_5=CreateSeuratObject(CountMatrix,project ="Normal Control CT5",min.cells=3,min.features=200)
srobj_5[["MTpercent"]]=PercentageFeatureSet(srobj_5,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/Normal Control/GSM4288845_CT6-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/Normal Control/GSM4288845_CT6-CD3-genes-barcodes-matrix")


srobj_6=CreateSeuratObject(CountMatrix,project ="Normal Control CT6",min.cells=3,min.features=200)
srobj_6[["MTpercent"]]=PercentageFeatureSet(srobj_6,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/Normal Control/GSM4288846_CT7-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/Normal Control/GSM4288846_CT7-CD3-genes-barcodes-matrix")


srobj_7=CreateSeuratObject(CountMatrix,project ="Normal Control CT7",min.cells=3,min.features=200)
srobj_7[["MTpercent"]]=PercentageFeatureSet(srobj_7,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/Normal Control/GSM4288847_CT8-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/Normal Control/GSM4288847_CT8-CD3-genes-barcodes-matrix")


srobj_8=CreateSeuratObject(CountMatrix,project ="Normal Control CT8",min.cells=3,min.features=200)
srobj_8[["MTpercent"]]=PercentageFeatureSet(srobj_8,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/GSE144469/Data/CD3/+CPI no colitis/GSM4288834_NC1-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/GSE144469/Data/CD3/+CPI no colitis/GSM4288834_NC1-CD3-genes-barcodes-matrix")


srobj_9=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC1",min.cells=3,min.features=200)
srobj_9[["MTpercent"]]=PercentageFeatureSet(srobj_9,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI no colitis/GSM4288834_NC1-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI no colitis/GSM4288834_NC1-CD3-genes-barcodes-matrix")


srobj_10=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC2",min.cells=3,min.features=200)
srobj_10[["MTpercent"]]=PercentageFeatureSet(srobj_10,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI no colitis/GSM4288836_NC3-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI no colitis/GSM4288836_NC3-CD3-genes-barcodes-matrix")


srobj_11=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC3",min.cells=3,min.features=200)
srobj_11[["MTpercent"]]=PercentageFeatureSet(srobj_11,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI no colitis/GSM4288837_NC4-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI no colitis/GSM4288837_NC4-CD3-genes-barcodes-matrix")


srobj_12=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC4",min.cells=3,min.features=200)
srobj_12[["MTpercent"]]=PercentageFeatureSet(srobj_12,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI no colitis/GSM4288838_NC5-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI no colitis/GSM4288838_NC5-CD3-genes-barcodes-matrix")


srobj_13=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC5",min.cells=3,min.features=200)
srobj_13[["MTpercent"]]=PercentageFeatureSet(srobj_13,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI no colitis/GSM4288839_NC6-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI no colitis/GSM4288839_NC6-CD3-genes-barcodes-matrix")


srobj_14=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC6",min.cells=3,min.features=200)
srobj_14[["MTpercent"]]=PercentageFeatureSet(srobj_14,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288826_C1-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288826_C1-CD3-genes-barcodes-matrix")


srobj_15=CreateSeuratObject(CountMatrix,project ="+CPI colitis C1",min.cells=3,min.features=200)
srobj_15[["MTpercent"]]=PercentageFeatureSet(srobj_15,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288827_C2-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288827_C2-CD3-genes-barcodes-matrix")


srobj_16=CreateSeuratObject(CountMatrix,project ="+CPI colitis C2",min.cells=3,min.features=200)
srobj_16[["MTpercent"]]=PercentageFeatureSet(srobj_16,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288828_C3-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288828_C3-CD3-genes-barcodes-matrix")


srobj_17=CreateSeuratObject(CountMatrix,project ="+CPI colitis C3",min.cells=3,min.features=200)
srobj_17[["MTpercent"]]=PercentageFeatureSet(srobj_17,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288829_C4-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288829_C4-CD3-genes-barcodes-matrix")


srobj_18=CreateSeuratObject(CountMatrix,project ="+CPI colitis C4",min.cells=3,min.features=200)
srobj_18[["MTpercent"]]=PercentageFeatureSet(srobj_18,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288830_C5-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288830_C5-CD3-genes-barcodes-matrix")


srobj_19=CreateSeuratObject(CountMatrix,project ="+CPI colitis C5",min.cells=3,min.features=200)
srobj_19[["MTpercent"]]=PercentageFeatureSet(srobj_19,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288831_C6-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288831_C6-CD3-genes-barcodes-matrix")


srobj_20=CreateSeuratObject(CountMatrix,project ="+CPI colitis C6",min.cells=3,min.features=200)
srobj_20[["MTpercent"]]=PercentageFeatureSet(srobj_20,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288832_C7-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288832_C7-CD3-genes-barcodes-matrix")


srobj_21=CreateSeuratObject(CountMatrix,project ="+CPI colitis C7",min.cells=3,min.features=200)
srobj_21[["MTpercent"]]=PercentageFeatureSet(srobj_21,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288833_C8-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288833_C8-CD3-genes-barcodes-matrix")


srobj_22=CreateSeuratObject(CountMatrix,project ="+CPI colitis C8",min.cells=3,min.features=200)
srobj_22[["MTpercent"]]=PercentageFeatureSet(srobj_22,pattern = "^MT-")


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
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/integrated_obj Data")
merged_obj <- merge(srobj_1, y = list(srobj_2, srobj_3, srobj_4, 
                                      srobj_5, srobj_6, srobj_7, 
                                      srobj_8, srobj_9,srobj_10, 
                                      srobj_11, srobj_12,srobj_13,
                                      srobj_14, srobj_15,srobj_16,
                                      srobj_17, srobj_18,srobj_19,
                                      srobj_20, srobj_21,srobj_22))


merged_obj=NormalizeData(merged_obj,normalization.method = "LogNormalize",scale.factor = 10000)


merged_obj=FindVariableFeatures(merged_obj,selection.method = "vst",nfeatures = 2000)
merged_obj=subset(merged_obj, features = VariableFeatures(merged_obj))


merged_obj = ScaleData(merged_obj,features = rownames(merged_obj))


merged_obj = RunPCA(merged_obj)


merged_obj <- IntegrateLayers(object = merged_obj,
                        method = CCAIntegration,
                        orig.reduction = "pca", 
                        new.reduction = "integrated.cca",
                        verbose = FALSE)


merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])


png(filename = "ElbowPlot.png",width = 10000,height=4000,units ="px",res = 600)
ElbowPlot(merged_obj,ndims = 50)
dev.off()


merged_obj1 = merged_obj


merged_obj1=FindNeighbors(merged_obj1,dims = 1:20)
merged_obj1=FindClusters(merged_obj1,resolution = 0.2)


merged_obj1=RunUMAP(merged_obj1,dims = 1:20)
png(filename = "UMAP1.png",width = 10000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1,reduction = "umap",label = TRUE)
dev.off()


merged_obj1 <- subset(merged_obj1, subset = seurat_clusters %in% c(10, 11), invert = TRUE)
png(filename = "UMAP2.png",width = 10000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1,reduction = "umap",label = TRUE)
dev.off()


umap_data <- Embeddings(merged_obj1, reduction = "umap")
cells_to_keep <- rownames(umap_data[umap_data[, "umap_1"] <= 10, ])
merged_obj1 <- subset(merged_obj1, cells = cells_to_keep)
png(filename = "UMAP3.png",width = 10000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1,reduction = "umap",label = TRUE)
dev.off()


merged_clusters <- as.character(Idents(merged_obj1))
merged_clusters[merged_clusters == "1" | merged_clusters == "8"] <- "1"
Idents(merged_obj1) <- merged_clusters
png(filename = "UMAP4.png",width = 10000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1,reduction = "umap",label = TRUE)
dev.off()


cluster_ids <- as.character(Idents(merged_obj1))
cluster_ids[cluster_ids == "9"] <- "8"
Idents(merged_obj1) <- cluster_ids
png(filename = "UMAP5.png",width = 10000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1,reduction = "umap",label = TRUE)
dev.off()


# UMAP projection of all integrated samples with cells colored by their Seurat clusters. Labels indicate cluster IDs.
desired_order <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")
Idents(merged_obj1) <- factor(Idents(merged_obj1), levels = desired_order)
png(filename = "UMAP_Clusters_AllSamples.png",width = 10000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1,reduction = "umap",label = TRUE)
dev.off()


# UMAP projection of all cells colored by their original sample ID (orig.ident).
png(filename = "UMAP_BySampleID.png",width = 10000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1, reduction = "umap", group.by = "orig.ident")
dev.off()


# UMAP projection split by original sample ID. Each panel shows the cells of one sample with their spatial distribution in the integrated UMAP space.
png(filename = "UMAP_SplitBySample.png",width = 25000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1, reduction = "umap", split.by = "orig.ident")
dev.off()
################################################################################ End integration and UMAP

merged_obj2 = merged_obj1

################################################################################ Start Find Marker
merged_obj2 <- JoinLayers(merged_obj2)
markers = FindAllMarkers(merged_obj2,min.pct = 0.1 , logfc.threshold = 0.1)
write.csv(markers,file="AllMarkers.csv")
#marker_2And3 = FindMarkers(merged_obj2,min.pct = 0.1 , logfc.threshold = 0.1, ident.1 = "2",ident.2 = "3")
################################################################################ End Find Marker

merged_obj3 = merged_obj2

################################################################################ Start annotation
# Dot Plot of Key Marker Genes Across All Clusters
png(filename = "DotPlot_Key_Markers_Across_All_Clusters.png",width = 10000,height=4000,units ="px",res = 600 )
DotPlot(merged_obj3, features = c("CD4", "CD40LG", "CD8A", "CD8B", "SELL", "CCR7", "IL7R", "CTLA4", 
                                  "FOXP3", "TNFRSF4", "TNFRSF18", "TIGIT", "IL2RA", "ICOS", "IFNG", 
                                  "STAT1", "STAT4", "TBX21", "GATA3", "STAT5", "STAT6", "IL17A", "RORA", 
                                  "LTB", "CD4A", "CD4B", "GNLY", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "NKG7")) + 
coord_flip()
dev.off()


important_markers <- c("CD4", "CD40LG", "CD8A", "CD8B", "SELL", "CCR7", "IL7R", "CTLA4", 
                       "FOXP3", "TNFRSF4", "TNFRSF18", "TIGIT", "IL2RA", "ICOS", "IFNG", 
                       "STAT1", "STAT4", "TBX21", "GATA3", "STAT5", "STAT6", "IL17A", "RORA", 
                       "LTB", "CD4A", "CD4B", "GNLY", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "NKG7")


# Filtered Key Markers from All Identified Cluster Markers
filtered_markers <- markers[markers$gene %in% important_markers, ]
write.csv(filtered_markers,file="Filtered_Markers.csv")
################################################################################ End annotation

