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
################################################################################ End Normal Control CT1 #2



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


srobj_3=subset(srobj_3,subset=nFeature_RNA>800 & nFeature_RNA<5000 & MTpercent<20 & nCount_RNA> 1800 & nCount_RNA<20000)
################################################################################ End Normal Control CT3 #3


################################################################################ Start +CPI no colitis NC1 #4
setwd("C:/Esmaeil/scRNA-seq/GSE144469/Data/CD3/+CPI no colitis/GSM4288834_NC1-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/GSE144469/Data/CD3/+CPI no colitis/GSM4288834_NC1-CD3-genes-barcodes-matrix")


srobj_4=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC1",min.cells=3,min.features=200)
srobj_4[["MTpercent"]]=PercentageFeatureSet(srobj_4,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_4@meta.data, aes(x = srobj_4$nCount_RNA, y = srobj_4$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 15, color = 'red') + 
  geom_vline(xintercept = 2000, color = 'red') +
  geom_vline(xintercept = 21000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_4@meta.data, aes(x=srobj_4$nFeature_RNA, y=srobj_4$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=15, color='red') + 
  geom_vline(xintercept=900, color='red') +
  geom_vline(xintercept=6000, color='red') 
dev.off()


srobj_4=subset(srobj_4,subset=nFeature_RNA>900 & nFeature_RNA<6000 & MTpercent<15 & nCount_RNA>2000 & nCount_RNA<21000)
################################################################################ End +CPI no colitis NC1 #4



################################################################################ Start +CPI no colitis NC2 #5
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI no colitis/GSM4288834_NC1-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI no colitis/GSM4288834_NC1-CD3-genes-barcodes-matrix")


srobj_5=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC2",min.cells=3,min.features=200)
srobj_5[["MTpercent"]]=PercentageFeatureSet(srobj_5,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_5@meta.data, aes(x = srobj_5$nCount_RNA, y = srobj_5$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 19, color = 'red') + 
  geom_vline(xintercept = 1600, color = 'red') +
  geom_vline(xintercept = 18000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_5@meta.data, aes(x=srobj_5$nFeature_RNA, y=srobj_5$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=19, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=4500, color='red') 
dev.off()


srobj_5=subset(srobj_5,subset=nFeature_RNA>800 & nFeature_RNA<4500 & MTpercent<19 & nCount_RNA>1600 & nCount_RNA<18000)
################################################################################ End +CPI no colitis NC2 #5



################################################################################ Start +CPI no colitis NC3 #6
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI no colitis/GSM4288836_NC3-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI no colitis/GSM4288836_NC3-CD3-genes-barcodes-matrix")


srobj_6=CreateSeuratObject(CountMatrix,project ="+CPI no colitis NC3",min.cells=3,min.features=200)
srobj_6[["MTpercent"]]=PercentageFeatureSet(srobj_6,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_6@meta.data, aes(x = srobj_6$nCount_RNA, y = srobj_6$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1600, color = 'red') +
  geom_vline(xintercept = 18000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_6@meta.data, aes(x=srobj_6$nFeature_RNA, y=srobj_6$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=4200, color='red') 
dev.off()


srobj_6=subset(srobj_6,subset=nFeature_RNA>800 & nFeature_RNA<4200 & MTpercent<20 & nCount_RNA>1600 & nCount_RNA<18000)
################################################################################ End +CPI no colitis NC3 #6



################################################################################ Start +CPI colitis C1 #7
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288826_C1-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288826_C1-CD3-genes-barcodes-matrix")


srobj_7=CreateSeuratObject(CountMatrix,project ="+CPI colitis C1",min.cells=3,min.features=200)
srobj_7[["MTpercent"]]=PercentageFeatureSet(srobj_7,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_7@meta.data, aes(x = srobj_7$nCount_RNA, y = srobj_7$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 25, color = 'red') + 
  geom_vline(xintercept = 1600, color = 'red') +
  geom_vline(xintercept = 30000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_7@meta.data, aes(x=srobj_7$nFeature_RNA, y=srobj_7$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=25, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=6200, color='red') 
dev.off()


srobj_7=subset(srobj_7,subset=nFeature_RNA>800 & nFeature_RNA<6200 & MTpercent<25 & nCount_RNA>1600 & nCount_RNA<30000)
################################################################################ End +CPI colitis C1 #7



################################################################################ Start +CPI colitis C2 #8
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288827_C2-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288827_C2-CD3-genes-barcodes-matrix")


srobj_8=CreateSeuratObject(CountMatrix,project ="+CPI colitis C2",min.cells=3,min.features=200)
srobj_8[["MTpercent"]]=PercentageFeatureSet(srobj_8,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_8@meta.data, aes(x = srobj_8$nCount_RNA, y = srobj_8$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1600, color = 'red') +
  geom_vline(xintercept = 37000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_8@meta.data, aes(x=srobj_8$nFeature_RNA, y=srobj_8$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=900, color='red') +
  geom_vline(xintercept=6200, color='red') 
dev.off()


srobj_8=subset(srobj_8,subset=nFeature_RNA>900 & nFeature_RNA<6200 & MTpercent<20 & nCount_RNA>1600 & nCount_RNA<37000)
################################################################################ End +CPI colitis C2 #8



################################################################################ Start +CPI colitis C3 #9
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288828_C3-CD3-genes-barcodes-matrix")
CountMatrix=Read10X("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Data/CD3/+CPI colitis/GSM4288828_C3-CD3-genes-barcodes-matrix")


srobj_9=CreateSeuratObject(CountMatrix,project ="+CPI colitis C3",min.cells=3,min.features=200)
srobj_9[["MTpercent"]]=PercentageFeatureSet(srobj_9,pattern = "^MT-")


png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj_9@meta.data, aes(x = srobj_9$nCount_RNA, y = srobj_9$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 20, color = 'red') + 
  geom_vline(xintercept = 1400, color = 'red') +
  geom_vline(xintercept = 37000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj_9@meta.data, aes(x=srobj_9$nFeature_RNA, y=srobj_9$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=20, color='red') + 
  geom_vline(xintercept=800, color='red') +
  geom_vline(xintercept=6200, color='red') 
dev.off()


srobj_9=subset(srobj_9,subset=nFeature_RNA>800 & nFeature_RNA<6200 & MTpercent<20 & nCount_RNA>1400 & nCount_RNA<37000)
################################################################################ End +CPI colitis C3 #9
