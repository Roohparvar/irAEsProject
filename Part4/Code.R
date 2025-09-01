################################################################################ Start library
library(Seurat)
library(SeuratObject)
library(ggplot2)
################################################################################ End library



################################################################################ Start P1_T #1
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part4/P1_T")

srobj=CreateSeuratObject(CountMatrix,project ="P1_T_Blood",min.cells=3,min.features=200)
srobj[["MTpercent"]]=PercentageFeatureSet(srobj,pattern = "^MT-")

setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part4/Quality Control/P1_T")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj@meta.data, aes(x = srobj$nCount_RNA, y = srobj$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 40, color = 'red') + 
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 22000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj@meta.data, aes(x=srobj$nFeature_RNA, y=srobj$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=40, color='red') + 
  geom_vline(xintercept=500, color='red') +
  geom_vline(xintercept=4000, color='red') 
dev.off()


srobj_P1_T=subset(srobj,subset=nFeature_RNA>500 & nFeature_RNA<4000 & MTpercent<40 & nCount_RNA>1000 & nCount_RNA<22000)
################################################################################ End P1_T #1



################################################################################ Start P1_TN #2
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part4/P1_TN")

rownames(CountMatrix) <- make.unique(rownames(CountMatrix))
colnames(CountMatrix) <- make.unique(colnames(CountMatrix))
rownames(CountMatrix) <- gsub("_", "-", rownames(CountMatrix))

srobj=CreateSeuratObject(CountMatrix,project ="P1_TN_Blood",min.cells=3,min.features=200)
srobj[["MTpercent"]]=PercentageFeatureSet(srobj,pattern = "^MT-")

setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part4/Quality Control/P1_TN")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj@meta.data, aes(x = srobj$nCount_RNA, y = srobj$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 1, color = 'red') + 
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 22000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj@meta.data, aes(x=srobj$nFeature_RNA, y=srobj$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=1, color='red') + 
  geom_vline(xintercept=500, color='red') +
  geom_vline(xintercept=4000, color='red') 
dev.off()


srobj_P1_TN=subset(srobj,subset=nFeature_RNA>500 & nFeature_RNA<4000 & MTpercent<40 & nCount_RNA>1000 & nCount_RNA<22000)
################################################################################ End P1_TN #2




################################################################################ Start P1_TN #2
CountMatrix=Read10X("C:/Esmaeil/irAEsProject/Backup/Part4/P6_T")

rownames(CountMatrix) <- make.unique(rownames(CountMatrix))
colnames(CountMatrix) <- make.unique(colnames(CountMatrix))
rownames(CountMatrix) <- gsub("_", "-", rownames(CountMatrix))

srobj=CreateSeuratObject(CountMatrix,project ="P6_T_Blood",min.cells=3,min.features=200)





install.packages("scCustomize")

# Load the package
library(scCustomize)

# Access the list of human mitochondrial Ensembl gene IDs
mt_genes <- ensembl_mito_id$Homo_sapiens_mito_ensembl

# Print the first few entries to verify
head(mt_genes)


srobj[["MTpercent"]] <- PercentageFeatureSet(srobj, features = mt_genes)

# Check the summary statistics
summary(srobj[["percent.mt"]])


# srobj[["MTpercent"]]=PercentageFeatureSet(srobj,pattern = "^MT-")

setwd("C:/Esmaeil/irAEsProject/irAEsProject/Part4/Quality Control/P6_T")

png(filename = "1.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data = srobj@meta.data, aes(x = srobj$nCount_RNA, y = srobj$MTpercent)) +
  geom_point(size = 2, color = 'blue') +
  labs(x = 'nCount_RNA', y = 'perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept = 1, color = 'red') + 
  geom_vline(xintercept = 1000, color = 'red') +
  geom_vline(xintercept = 22000, color = 'red')
dev.off()


png(filename = "2.png",width = 10000,height=4000,units ="px",res = 600)
ggplot(data=srobj@meta.data, aes(x=srobj$nFeature_RNA, y=srobj$MTpercent)) +
  geom_point(size=2, color = 'blue') +
  labs(x='nFeature_RNA', y='perc. mito') +
  scale_x_log10() +
  geom_hline(yintercept=1, color='red') + 
  geom_vline(xintercept=500, color='red') +
  geom_vline(xintercept=4000, color='red') 
dev.off()


srobj_P2_TN=subset(srobj,subset=nFeature_RNA>500 & nFeature_RNA<4000 & MTpercent<40 & nCount_RNA>1000 & nCount_RNA<22000)
################################################################################ End P1_TN #6

















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
png(filename = "UMAP1.png",width = 10000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1,label = TRUE)
dev.off()


# Remove clusters 11
merged_obj1 <- subset(merged_obj1, subset = seurat_clusters %in% c(11), invert = TRUE)
png(filename = "UMAP2.png",width = 10000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1,label = TRUE)
dev.off()


# UMAP projection of all cells colored by their original sample ID (orig.ident).
png(filename = "UMAP3.png",width = 7000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1, group.by = "orig.ident")
dev.off()


# UMAP projection split by original sample ID. Each panel shows the cells of one sample with their spatial distribution in the integrated UMAP space.
png(filename = "UMAP4.png",width = 28000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1, split.by = "orig.ident")
dev.off()


png(filename = "FeaturePlot.png", width = 10000, height = 4000, units = "px", res = 600)
FeaturePlot(merged_obj1, features = c("CD3D", "CD3E", "CD3G", "LYZ", "CD79A", "CD19"))
dev.off()


# Remove clusters 10
merged_obj1 <- subset(merged_obj1, subset = seurat_clusters %in% c(10), invert = TRUE)
png(filename = "UMAP5.png",width = 10000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1,label = TRUE)
dev.off()



# Remove cells with an x-axis value less than -11 and y-axis value less than -1.
umap_coord <- merged_obj1@reductions[["umap"]]@cell.embeddings
cells_to_keep <- rownames(umap_coord[!(umap_coord[, "umap_1"] < -11 & umap_coord[, "umap_2"] < -1), ])
merged_obj1 <- subset(merged_obj1, cells = cells_to_keep)
png(filename = "UMAP6.png",width = 10000,height=4000,units ="px",res = 600 )
DimPlot(merged_obj1, label = TRUE)
dev.off()
################################################################################ End UMAP

merged_obj2 = merged_obj1

# 4
# saveRDS(file = "merged_obj2",merged_obj2)
# The Seurat object obtained after UMAP


################################################################################ Start Extracting and saving Seurat objects for each sample

samples <- unique(merged_obj2$orig.ident)

for (i in seq_along(samples)) {
  sample_name <- samples[i]
  seurat_obj <- subset(merged_obj2, subset = orig.ident == sample_name)
  assign(paste0("srobj_", i), seurat_obj)
}


# 5
# 5_GSE144469_seurat_objs
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

