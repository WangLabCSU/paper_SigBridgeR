library(Seurat)
library(ggplot2)
library(patchwork)
library(zeallot)
library(dplyr)

setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_binary/brca/TNBC")

data_path = "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_data/brca"


GSE162228_pheno = qs::qread(file.path(data_path, "brca_pheno_GSE162228.qs"))

seurat_merged = qs::qread("GSE162228_tnbc_merged_seurat.qs")

set.seed(123)
my_color = randomcoloR::distinctColorPalette(
  length(unique(seurat_merged$seurat_clusters)),
  runTsne = TRUE
)

# 更新Seurat对象
# counts <- LayerData(seurat_merged, assay = "RNA", layer = "counts")
# counts <- GetAssayData(seurat_merged, assay = "RNA", slot = "counts")

# metadata <- seurat_merged@meta.data

# seurat_merged = CreateSeuratObject(
#     counts = counts,
#     meta.data = metadata
# ) %>%
#     NormalizeData() %>%
#     FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#     ScaleData() %>%
#     RunPCA(features = VariableFeatures(.)) %>%
#     RunUMAP(dims = 1:30)

seurat_tnbc_tumor = readRDS(
  "/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_TNBCTum.rds"
)
seurat_merged$is_tumor = ifelse(
  colnames(seurat_merged) %in% rownames(seurat_tnbc_tumor@meta.data),
  "TRUE",
  "FALSE"
)

umap_cluster = DimPlot(
  seurat_merged,
  reduction = "umap",
  group.by = "seurat_clusters",
  pt.size = 0.05,
  label = FALSE,
  label.size = 5,
  cols = my_color
)

umap_tumor = DimPlot(
  seurat_merged,
  reduction = "umap",
  group.by = "is_tumor",
  pt.size = 0.05,
  cols = c("FALSE" = "#386c9b", "TRUE" = "#ff3333")
)


c(
  umap_scissor,
  umap_scpas,
  umap_scpp,
  umap_scab,
  umap_lp_sgl,
  umap_degas
) %<-%
  purrr::map(
    c(
      "scissor",
      "scPAS",
      "scPP",
      'scAB',
      "LP_SGL",
      "DEGAS"
    ),
    ~ Seurat::DimPlot(
      seurat_merged,
      group.by = .x,
      label.size = 8,
      pt.size = 0.05,
      reduction = "umap",
      cols = c(
        "Other" = "#CECECE",
        "Neutral" = "#CECECE",
        "Positive" = "#ff3333",
        "Negative" = "#386c9b"
      )
    ) +
      ggplot2::ggtitle(.x)
  )

umaps =
  umap_cluster +
  umap_tumor +
  umap_scab +
  umap_scissor +
  umap_scpas +
  umap_scpp +
  umap_degas +
  umap_lp_sgl +
  plot_layout(ncol = 3)

ggsave(
  umaps,
  filename = "GSE162228_UMAP_screened.png",
  width = 15,
  height = 14,
  dpi = 400
)
