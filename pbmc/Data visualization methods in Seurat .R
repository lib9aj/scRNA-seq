
###Data visualization methods in Seurat

#We’ll demonstrate visualization techniques in Seurat using our previously computed Seurat object from the 2,700 PBMC tutorial. 
#You can download this dataset from SeuratData

SeuratData::InstallData("pbmc3k")

library (Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)

pbmc3k.final <- LoadData("pbmc3k", type = "pbmc3k.final")
pbmc3k.final$groups <- sample(c ("group1", "group2"), size = ncol(pbmc3k.final), replace = TRUE)
features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")
pbmc3k.final
DimPlot(pbmc3k.final, reduction = "umap")

###Five visualizations of marker feature expression

# 1.Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(pbmc3k.final, features = features, ncol = 2)

# 2.Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(pbmc3k.final, features = features)

# 3.Feature plot - visualize feature expression in low-dimensional space
FeaturePlot(pbmc3k.final, features = features)

# 4.Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(pbmc3k.final, features = features) + RotatedAxis()

# 5.Single cell heatmap of feature expression
DoHeatmap(subset(pbmc3k.final, downsample = 100), features = features, size = 3)



