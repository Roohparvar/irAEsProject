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