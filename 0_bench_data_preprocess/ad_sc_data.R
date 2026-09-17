ad_sc <- readRDS(
  "/home/data/sigbridger/benchmark_data/ad/Nat_Neurosci_10.1038_s41593-019-0539-4/adsn_seuV5.rds"
)

# 已有的预处理：NormalizeData / FindVariableFeatures / ScaleData / RunPCA(30 dims) / RunUMAP
# 补充缺失的步骤：邻居图与聚类
ad_sc <- Seurat::FindNeighbors(ad_sc, dims = 1:30)
ad_sc <- Seurat::FindClusters(ad_sc, resolution = 0.8)

table(ad_sc$cellType)

ad_sc <- ad_sc[, ad_sc$cellType != "doublet"]

table(ad_sc$cellType)


qs::qsave(
  ad_sc,
  "/home/data/sigbridger/benchmark_data/ad/ad_GSE138852_seurat.qs",
  nthreads = 4L
)
