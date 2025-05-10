################################################################################ Start library
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(openxlsx)
################################################################################ End library



################################################################################ Start Load Seurat objects from Dataset 1: GSE144469
folder_path1 <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 3/Data/Part 1_GSE144469_Seurat_Objs"

rds_files1 <- list.files(folder_path1, pattern = "\\.rds$", full.names = TRUE)

for (i in seq_along(rds_files1)) {
  obj <- readRDS(rds_files1[i])
  assign(paste0("srob", i), obj)
}
################################################################################ End Load Seurat objects from Dataset 1: GSE144469



################################################################################ Start Load Seurat objects from Dataset 2: GSE206299
folder_path2 <- "C:/Esmaeil/scRNA-seq/Backup of Local Data and Files Not on GitHub/Part 3/Data/Part 2_GSE206299_Seurat_Objs"

rds_files2 <- list.files(folder_path2, pattern = "\\.rds$", full.names = TRUE)

offset <- length(rds_files1)  
for (i in seq_along(rds_files2)) {
  obj <- readRDS(rds_files2[i])
  assign(paste0("srob", i + offset), obj)
}
################################################################################ End Load Seurat objects from Dataset 2: GSE206299



################################################################################ Start integration
setwd("C:/Esmaeil/scRNA-seq/Single-Cell-Pipeline-in-R/Part3")

merged_obj <- merge(srob1, y = list(srob2, srob3, srob4, 
                                    srob5, srob6, srob7, 
                                    srob8, srob9,srob10, 
                                    srob11, srob12,srob13,
                                    srob14, srob15,srob16,
                                    srob17, srob18,srob19,
                                    srob20, srob21,srob22,
                                    srob23, srob24,srob25,
                                    srob26, srob27,srob28,
                                    srob29, srob30,srob31,
                                    srob32, srob32,srob33,
                                    srob34, srob35,srob36,
                                    srob37, srob38,srob39,
                                    srob40, srob41,srob42,
                                    srob43, srob44,srob45,
                                    srob46, srob47,srob48,
                                    srob49, srob50,srob51,srob52))

gc()
merged_obj=NormalizeData(merged_obj,normalization.method = "LogNormalize",scale.factor = 10000)

gc()

merged_obj=FindVariableFeatures(merged_obj,selection.method = "vst",nfeatures = 2000)
#merged_obj=subset(merged_obj, features = VariableFeatures(merged_obj))
gc()

merged_obj = ScaleData(merged_obj)
gc()

merged_obj = RunPCA(merged_obj)
gc()

# 1
# saveRDS(file = "merged_obj",merged_obj)
# The Seurat object obtained after RunPCA and before IntegrateLayers
gc()
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

merged_obj1@meta.data <- merged_obj1@meta.data %>% mutate(dataset = case_when(
  grepl(x = merged_obj1$orig.ident, pattern = "CPI") ~ "GSE144469",
  grepl(x = merged_obj1$orig.ident, pattern = "Normal") ~ "GSE144469",
  grepl(x = merged_obj1$orig.ident, pattern = "irColitis") ~ "GSE206299",
  grepl(x = merged_obj1$orig.ident, pattern = "Healthy") ~ "GSE206299",
  grepl(x = merged_obj1$orig.ident, pattern = "Therapy") ~ "GSE206299"
))


merged_obj1=FindNeighbors(merged_obj1,dims = 1:30,reduction = "integrated.cca")
merged_obj1=FindClusters(merged_obj1,resolution = 0.17)


merged_obj1=RunUMAP(merged_obj1,dims = 1:30, reduction = "integrated.cca")
png(filename = "UMAP1.png", width = 8000, height = 4000, units = "px", res = 1200)
DimPlot(merged_obj1, label = TRUE)
dev.off()


png(filename = "UMAP2.png",width = 20000,height=9000,units ="px",res = 600 )
DimPlot(merged_obj1, split.by = "dataset") + 
  theme(
    strip.text = element_text(size = 35),   # اندازه عنوان پنل ها
    axis.text = element_text(size = 20),    # اندازه اعداد محور
    axis.title = element_text(size = 20),   # عنوان محورها
    legend.text = element_text(size = 17) 
  ) +
  guides(color = guide_legend(override.aes = list(size = 8))) 
dev.off()

################################################################################ End UMAP

merged_obj2 = merged_obj1

# 3
# saveRDS(file = "merged_obj2",merged_obj2)
# The Seurat object obtained after UMAP

################################################################################ Start Find Marker
merged_obj2 <- JoinLayers(merged_obj2)
markers = FindAllMarkers(merged_obj2,min.pct = 0.1 , logfc.threshold = 0.1)
write.csv(markers,file="AllMarkers.csv")
################################################################################ ٍدی Find Marker

# 4
# saveRDS(file = "merged_obj2",merged_obj2)
# The Seurat object obtained after  FindAllMarkers

################################################################################ Start finding Significant Genes and saving each cluster's into a separate sheet
df <- read.csv("AllMarkers.csv")
filtered_df <- df %>% filter(
  p_val < 0.05,
  pct.1 > 0.1,
  (pct.1 - pct.2) > 0.1
)

wb <- createWorkbook()

clusters <- unique(filtered_df$cluster)

for (i in seq_along(clusters)) {
  cluster_value <- clusters[i]
  
  cluster_data <- filtered_df %>% filter(cluster == cluster_value)
  
  sheet_name <- paste0("Cluster_", cluster_value)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, cluster_data)
}

saveWorkbook(wb, "Significant_Genes.xlsx", overwrite = TRUE)
################################################################################ End finding Significant Genes and saving each cluster's into a separate sheet



################################################################################ Start Dot Plot of Key Marker Genes
png(filename = "DotPlot.png",width = 10000,height=4000,units ="px",res = 600 )
DotPlot(merged_obj2, features = c("CD4", "CD40LG", "CD8A", "CD8B", "SELL", "CCR7", "IL7R", "CTLA4", 
                                  "FOXP3", "TNFRSF4", "TNFRSF18", "TIGIT", "IL2RA", "ICOS", "IFNG", 
                                  "STAT1", "STAT4", "TBX21", "GATA3", "STAT5", "STAT6", "IL17A", "RORA", 
                                  "LTB", "CD4A", "CD4B", "GNLY", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "NKG7")) + 
coord_flip()
dev.off()
################################################################################ End Dot Plot of Key Marker Genes