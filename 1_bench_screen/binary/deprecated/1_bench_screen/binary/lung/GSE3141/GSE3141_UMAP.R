library(Seurat)
library(ggplot2)
library(patchwork)
library(zeallot)
library(dplyr)

setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_data/lung/GSE3141")

set.seed(123)

seurat_merged = qs::qread("GSE3141_luad_merged_seurat.qs", nthreads = 4)
my_color = randomcoloR::distinctColorPalette(
  k = length(unique(seurat_merged$seurat_clusters)),
  runTsne = T
)


c(umap_cnv, umap_tissue) %<-%
  purrr::map(
    c("cnv_status", "Tissue"),
    ~ Seurat::DimPlot(
      seurat_merged,
      reduction = "umap",
      group.by = .x,
      pt.size = 0.05,
      cols = c(
        "tumor" = "#ff3333",
        "Tumor" = "#ff3333",
        "normal" = "#386c9b",
        "Normal" = "#386c9b"
      )
    ) +
      ggplot2::ggtitle(.x)
  )

c(
  umap_scissor,
  umap_scpas,
  # umap_scpp,
  umap_scab
) %<-%
  purrr::map(
    c(
      "scissor",
      "scPAS",
      # "scPP",
      "scAB"
    ),
    ~ Seurat::DimPlot(
      seurat_merged,
      group.by = .x,
      reduction = "umap",
      pt.size = 0.05,
      cols = c(
        "Other" = "#CECECE",
        "Neutral" = "#CECECE",
        "Positive" = "#ff3333",
        "Negative" = "#386c9b"
      )
    ) +
      ggplot2::ggtitle(.x)
  )

c(umap_cluster, umap_sample, umap_celltype) %<-%
  purrr::map(
    c(
      "seurat_clusters",
      "Sample",
      "Celltype"
    ),
    ~ Seurat::DimPlot(
      seurat_merged,
      group.by = .x,
      pt.size = 0.05,
      reduction = "umap",
      cols = my_color
    ) +
      ggplot2::ggtitle(.x)
  )

umaps = umap_cluster +
  umap_sample +
  umap_celltype +
  umap_scissor +
  umap_scpas +
  umap_scab +
  umap_cnv +
  umap_tissue +
  plot_layout(ncol = 3)

ggplot2::ggsave(
  umaps,
  filename = "GSE3141_UMAP_screened.png",
  width = 14,
  height = 10,
  dpi = 400
)
