
library(Seurat)
library(dplyr)
library(hdf5r)
library(jsonlite)
library(harmony)


# Define .h5 file paths
h5_paths <- list(
  "dpc2.5" = "/Users/lib9aj/git/scRNA-seq/dpc2.5_3.5_4_4.5_Haixiang_Sun/dpc2.5/dpc2.5_filtered_feature_bc_matrix.h5",
  "dpc3.5" = "/Users/lib9aj/git/scRNA-seq/dpc2.5_3.5_4_4.5_Haixiang_Sun/dpc3.5/dpc3.5_filtered_feature_bc_matrix.h5",
  "dpc4" = "/Users/lib9aj/git/scRNA-seq/dpc2.5_3.5_4_4.5_Haixiang_Sun/dpc4/dpc4_filtered_feature_bc_matrix.h5",
  "dpc4.5" = "/Users/lib9aj/git/scRNA-seq/dpc2.5_3.5_4_4.5_Haixiang_Sun/dpc4.5/dpc4.5_filtered_feature_bc_matrix.h5"
)


# 读取每个时间点的数据并创建 Seurat 对象
seurat_list <- lapply(names(h5_paths), function(day) {
  data <- Read10X_h5(h5_paths[[day]])
  seurat_obj <- CreateSeuratObject(counts = data, project = day)
  seurat_obj$day <- day  # 添加时间点元数据
  
  # 计算线粒体基因百分比
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  
  # 质量控制，去除低质量细胞
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & percent.mt < 15)
  
  print(paste(day, "cell numbers:", ncol(seurat_obj)))  # 打印每个数据集的细胞数量
  
  return(seurat_obj)
})

# 合并 Seurat 对象
combined_seurat <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(h5_paths))

# 检查合并后每个时间点的细胞数量
table(combined_seurat$day)

# 添加 run_id 信息
cell_ids <- colnames(combined_seurat)
run_ids <- sapply(strsplit(cell_ids, "_"), `[`, 1)
combined_seurat$run_id <- run_ids

# 数据标准化
combined_seurat <- NormalizeData(combined_seurat)

# 寻找高变基因
combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 2000)

# 数据缩放
combined_seurat <- ScaleData(combined_seurat)

# 主成分分析（PCA）
combined_seurat <- RunPCA(combined_seurat, features = VariableFeatures(object = combined_seurat))

# 批次效应校正（Harmony）
combined_seurat <- RunHarmony(combined_seurat, "day")

# 提取合并后的元数据
meta_data <- combined_seurat@meta.data

# 统计合并后的详细信息
detailed_counts <- meta_data %>%
  group_by(day, run_id) %>%
  summarise(
    cell_count = n(),
    total_features = sum(nFeature_RNA, na.rm = TRUE),
    total_counts = sum(nCount_RNA, na.rm = TRUE)
  ) %>%
  group_by(day) %>%
  summarise(
    run_count = n(),
    run_ids = paste(run_id, collapse = ", "),
    total_cells = sum(cell_count),
    total_features = sum(total_features),
    avg_cells_per_run = mean(cell_count),
    avg_features_per_run = mean(total_features)
  )

# 打印详细统计信息
print(detailed_counts)

# UMAP降维
combined_seurat <- RunUMAP(combined_seurat, reduction = "harmony", dims = 1:10)

# 细胞聚类
combined_seurat <- FindNeighbors(combined_seurat, reduction = "harmony", dims = 1:10)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)

# 绘制 UMAP 图，按天数和 run ID 着色

DimPlot(combined_seurat, reduction = "umap", group.by = "day", raster = FALSE)
DimPlot(combined_seurat, reduction = "umap", group.by = "run_id", raster = FALSE)












