out_dir <- "/home/data/sigbridger"
data <- anndataR::read_h5ad(data_dir)


seurat <- Seurat::LoadXenium(data_dir, fov = "fov")
seurat <- subset(
  seurat,
  nCount_Xenium > 0L &
    nCount_Xenium < 50000L &
    nFeature_Xenium > 0L &
    nFeature_Xenium < 6000L
)

Seurat::DefaultAssay(seurat) # "Xenium"

xenium.obj <- Seurat::SCTransform(seurat, assay = "Xenium")
xenium.obj <- Seurat::RunPCA(
  xenium.obj,
  npcs = 30,
  features = rownames(xenium.obj)
)
xenium.obj <- Seurat::RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- Seurat::FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- Seurat::FindClusters(xenium.obj, resolution = 0.3)

qs::qsave(xenium.obj, file.path(out_dir, "xenium.obj.qs"), nthreads = 4L)

Seurat::VlnPlot(
  xenium.obj,
  features = c("nFeatures_Xenium", "nCount_Xenium"),
  col = 2,
  pt.size = 0
)

Seurat::ImageDimPlot(
  xenium.obj,
  molecules = c("Sst", "Gad1", "Pvalb", "Gfap"),
  nmols = 20000
)

xenium.obj <- qs::qread(file.path(out_dir, "xenium.obj.qs"), nthreads = 4L)

# ----------------------------------------------------------------------------------------------------
# 1. 为每个样本创建标准目录结构
samples <- c(
  "GSM8443449_PDAC-p1",
  "GSM8443450_PDAC-p2",
  "GSM8443451_PDAC-p3",
  "GSM8443452_PDAC-p4",
  "GSM8443453_PDAC-p5"
)

for (sample in samples) {
  # 创建标准目录
  sample_dir <- file.path("/home/data/sigbridger/spatial/GSE274103", sample)
  dir.create(sample_dir, showWarnings = FALSE)

  # 移动/复制H5矩阵
  h5_file <- file.path(
    "/home/data/sigbridger/spatial/GSE274103",
    paste0(sample, "_filtered_feature_bc_matrix.h5")
  )
  file.copy(h5_file, file.path(sample_dir, "filtered_feature_bc_matrix.h5"))

  # 解压spatial文件（.tar或.tar.gz）
  spatial_file <- file.path(
    "/home/data/sigbridger/spatial/GSE274103",
    paste0(sample, "_spatial.tar")
  )
  if (!file.exists(spatial_file)) {
    spatial_file <- paste0(spatial_file, ".gz") # 尝试.gz扩展名
  }
  untar(spatial_file, exdir = sample_dir) # 自动解压到sample_dir/spatial/
}

# rename directory as spatial
for (sample in samples) {
  sample_dir <- file.path("/home/data/sigbridger/spatial/GSE274103", sample)
  dirs <- list.dirs(sample_dir, full.names = TRUE, recursive = TRUE)

  file.rename(dirs[[2]], file.path(dirname(dirs[[2]]), "spatial"))
}

# R中执行
file.rename(
  from = "/home/data/sigbridger/spatial/GSE274103/GSM8443449_PDAC-p1/C30_spatial",
  to = "/home/data/sigbridger/spatial/GSE274103/GSM8443449_PDAC-p1/spatial"
)

# 读取数据
library(Seurat)
spatial_list <- lapply(samples, function(sample) {
  Load10X_Spatial(
    data.dir = file.path("/home/data/sigbridger/spatial/GSE274103", sample),
    slice = sample,
    assay = "Spatial"
  )
})

spatial_list <- lapply(spatial_list, function(seurat) {
  seurat <- subset(
    seurat,
    nCount_Spatial > 0L &
      nFeature_Spatial > 0
  )

  SCTransform(seurat, assay = "Spatial")
})

spatials <- merge(
  x = spatial_list[[1]],
  y = unlist(spatial_list[-1]),
  add.cell.ids = samples,
  collapse = FALSE,
  merge.data = TRUE,
  merge.dr = FALSE,
  project = "Integrated Spatial Seurat"
)

# DefaultAssay(spatials) <- "SCT"
VariableFeatures(spatials) <- lapply(spatial_list, function(seurat) {
  VariableFeatures(seurat)
}) %>%
  unlist() %>%
  unique()
spatials <- RunPCA(spatials, verbose = FALSE)
spatials <- FindNeighbors(spatials, dims = 1:30)
spatials <- FindClusters(spatials, verbose = FALSE)
spatials <- RunUMAP(spatials, dims = 1:30)

qs::qsave(spatials, "/home/data/sigbridger/spatials.qs", nthreads = 4L)

spatials <- qs::qread("/home/data/sigbridger/spatials.qs", nthreads = 4L)

SpatialDimPlot(spatials)

DefaultAssay(spatials)


bulk <- qs::qread(
  "/home/data/mutational_signature_sc_screen/bulk/TCGA-PAAD.exp.count.qs"
)

library(UCSCXenaShiny)
tcga_surv = load_data("tcga_surv")
# tcga_pheno = load_data("tcga_clinical")
# tcga_pheno_fine = load_data("tcga_clinical_fine")

# table(tcga_pheno$type)

#  ACC BLCA BRCA CESC CHOL COAD DLBC ESCA  GBM HNSC KICH KIRC KIRP LAML  LGG LIHC LUAD LUSC MESO   OV PAAD PCPG PRAD READ SARC SKCM STAD TGCT THCA
#   92  442 1236  312   45  554   48  204  604  604   91  952  352  200  533  439  661  624   87  604  196  187  569  183  271  479  511  139  580
# THYM UCEC  UCS  UVM
#  126  591   57   80

paad_samples <- colnames(bulk)

# ? survival phenotype
tcga_paad_surv = filter(tcga_surv, sample %in% paad_samples) %>%
  select(sample, OS.time, OS) %>%
  rename(time = OS.time, status = OS) %>%
  filter(status != "NA", time != "NA", !is.na(status), !is.na(time)) %>%
  tibble::column_to_rownames(var = "sample")

bulk <- bulk[, rownames(tcga_paad_surv)]

setThreads(8L)

# spatials@assays$RNA <- spatials@assays$SCT

res <- Screen(
  matched_bulk = bulk,
  sc_data = spatials,
  phenotype = tcga_paad_surv,
  label_type = NULL,
  phenotype_class = c("survival"),
  screen_method = "Scissor",
  assay = "SCT"
)

qs::qsave(res, "/home/data/sigbridger/spatials_res.qs", nthreads = 4L)

spatials_res <- qs::qread(
  "/home/data/sigbridger/spatials_res.qs",
  nthreads = 4L
)

Seurat::ImageDimPlot(spatials_res$scRNA_data)

positive_cell <- colnames(spatials_res$scRNA_data)[
  spatials_res$scRNA_data$scissor == "Positive"
]

p <- Seurat::SpatialDimPlot(
  spatials_res$scRNA_data,
  cells.highlight = positive_cell,
  ncol = 3L,
  # 降低背景切片透明度（0=完全透明，1=不透明），推荐 0.2~0.4
  image.alpha = 0.5,
  # 增强高亮细胞的视觉突出性
  cols.highlight = c("#ff3333", "#CECECE"), # 使用高对比度颜色
  alpha = c(0.8, 1) # c(非高亮点透明度, 高亮点透明度)
) &
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(
  "vignettes/example_figures/spatial_dim_plot.png",
  p,
  width = 7
)


seurat <- spatials_res$scRNA_data

data <- SeuratObject::LayerData(seurat)

adata <- anndataR::as_AnnData(
  x = seurat,
  x_mapping = "data",
  output_class = "ReticulateAnnData"
)

seurat_4 <- SCAnnotate(
  seurat,
  model = "Adult_Human_PancreaticIslet.pkl",
  method = "CellTypist",
  majority_voting = TRUE,
  download = TRUE,
  celltypist_tools = "/home/yyx/R/Project/R_code/SigBridgeR/inst/python/73-CellTypistAnnotate.py"
)

# Warning: No venv found in "~/.virtualenvs", "~/.venvs", "./venv", and "./.venv", return empty virtual environment result
# ℹ [2026/02/09 12:20:06] Existing environment "r-reticulate-celltypist" found
# ℹ [2026/02/09 12:20:06] [CellTypist] Start annotating cell types
# Warning: No venv found in , return empty virtual environment result
# ℹ [2026/02/09 12:20:09] Using Conda Env: "r-reticulate-celltypist"
# Warning: Skipping Layer "scale.data" with unexpected dimensions
# ℹ Expected [17716, 23436], got [723, 23436]
# ⚠️ Warning: invalid expression matrix, expect ALL genes and log1p normalized expression to 10000 counts per cell. The prediction result may not be accurate
# 🔬 Input data has 23436 cells and 17716 genes
# 🔗 Matching reference genes in the model
# 🧬 1617 features used for prediction
# ⚖️ Scaling input data
# 🖋️ Predicting labels
# ✅ Prediction done!
# 👀 Can not detect a neighborhood graph, will construct one before the over-clustering
# ⛓️ Over-clustering input data with resolution set to 15
# 🗳️ Majority voting the predictions
# ✅ Majority voting done!
# ℹ [2026/02/09 12:21:41] Annotation done

qs::qsave(
  seurat_4,
  "/home/data/sigbridger/celltypist_spatial_res.qs",
  nthreads = 4L
)

seurat_5 <- Seurat::PrepSCTFindMarkers(seurat)

seurat_5 <- SCAnnotate(
  seurat_5,
  #   method = "mLLMCelltype",
  models = c("deepseek-chat"),
  api_keys = list(
    deepseek = "sk-02cd71d0b3bd40479e5a6c41e84854d5"
  )
)

qs::qsave(
  seurat_5,
  "/home/data/sigbridger/mllmcelltype_spatial_res.qs",
  nthreads = 4L
)

spatials <- qs::qread("/home/data/sigbridger/mllmcelltype_spatial_res.qs")

# ℹ [2026/02/09 11:26:14] [mLLMCelltype] Start annotating cell types
# ℹ [2026/02/09 11:26:14] Find marker genes for each clusters
# Calculating cluster 0
# Calculating cluster 1
# Calculating cluster 2
# Calculating cluster 3
# Calculating cluster 4
# Calculating cluster 5
# Calculating cluster 6
# Calculating cluster 7
# Calculating cluster 8
# Calculating cluster 9
# Calculating cluster 10
# Calculating cluster 11
# Calculating cluster 12
# Calculating cluster 13
# Calculating cluster 14
# Calculating cluster 15
# Calculating cluster 16
# Calculating cluster 17
# Calculating cluster 18
# Calculating cluster 19
# Calculating cluster 20
# Calculating cluster 21
# Calculating cluster 22
# Calculating cluster 23
# Calculating cluster 24
# ℹ [2026/02/09 11:28:42] Large language models cell type Annotating
# ###
# # LLM Output
# ###
# ✔ [2026/02/09 11:28:46] Annotation Finished

merged <- MergeResult(seurat_5, seurat_4)

cell_type_names <- unique(seurat_5$cell_type2)
colors <- randomcoloR::distinctColorPalette(
  length(cell_type_names),
  runTsne = TRUE
)
names(colors) <- cell_type_names

seurat_5$cell_type2 <- seurat_5$scissor
seurat_5$cell_type2[
  seurat_5$scissor != "Positive"
] <- seurat_5$mllmcelltype_cell_type[seurat_5$scissor != "Positive"]

p <- Seurat::SpatialDimPlot(
  seurat_5,
  group.by = "cell_type2",
  ncol = 3L,
  cols = colors,
  # 降低背景切片透明度（0=完全透明，1=不透明），推荐 0.2~0.4
  image.alpha = 0.5,
  # 增强高亮细胞的视觉突出性
  cols.highlight = c("#ff3333", "#CECECE"), # 使用高对比度颜色
  alpha = c(0.8, 1), # c(非高亮点透明度, 高亮点透明度)
  stroke = 0
)

table(seurat_5$mllmcelltype_cell_type[seurat_5$scissor == "Positive"])

#     Acinar cells           Adipocytes              B cells         Chondrocytes    Endothelial cells
#              106                  100                   98                  566                   27
#      Enterocytes     Epithelial cells          Fibroblasts        Gastric cells         Goblet cells
#              106                  406                  319                   23                    9
#    Keratinocytes          Macrophages    Mesothelial cells Neuroendocrine cells  Smooth muscle cells
#                7                   17                    2                   43                   54

p2 <- Seurat::SpatialDimPlot(
  spatials,
  group.by = "mllmcelltype_cell_type",
  ncol = 3
)


cell_total <- table(spatials$mllmcelltype_cell_type)
cell_positive <- table(spatials$mllmcelltype_cell_type[
  spatials$scissor == "Positive"
])

ratios <- cell_positive[match(names(cell_total), names(cell_positive))] /
  cell_total
ratios[is.na(ratios)] <- 0

result <- sort(round(ratios, 4), decreasing = TRUE)
