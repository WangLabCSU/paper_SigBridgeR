library(Seurat)
library(ggplot2)
library(patchwork)
library(zeallot)
library(dplyr)

# setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_data//ov/GSE165897")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

seurat_merged = qs::qread("GSE9891_hgsoc_merged_seurat.qs", nthreads = 4)
set.seed(123)
my_color = randomcoloR::distinctColorPalette(
  k = 27,
  runTsne = TRUE
)


umap_treatment = DimPlot(
  seurat_merged,
  reduction = "umap",
  group.by = "treatment_phase",
  pt.size = 0.5,
  label.size = 7,
  cols = c(
    "#ff3333",
    "#386c9b"
  )
)

c(
  umap_scissor,
  umap_scpas,
  umap_scpp,
  umap_scab,
  umap_degas,
  umap_lp_sgl,
  umap_pipet
) %<-%
  purrr::map(
    c("scissor", "scPAS", "scPP", "scAB", "DEGAS", "LP_SGL", "PIPET"),
    ~ Seurat::DimPlot(
      seurat_merged,
      group.by = .x,
      pt.size = 0.5,
      reduction = "umap",
      cols = c(
        "Neutral" = "#CECECE",
        "Positive" = "#ff3333",
        "Negative" = "#386c9b",
        "Other" = "#CECECE"
      )
    ) +
      ggplot2::ggtitle(.x)
  )

c(umap_cluster, umap_location, umap_celltype, umap_subtype) %<-%
  purrr::map(
    c(
      "seurat_clusters",
      "anatomical_location",
      "cell_type",
      "cell_subtype"
    ),
    ~ Seurat::DimPlot(
      seurat_merged,
      group.by = .x,
      pt.size = 0.5,
      reduction = "umap",
      cols = my_color
    ) +
      ggplot2::ggtitle(.x)
  )

umaps = umap_cluster +
  umap_location +
  umap_treatment +
  umap_scissor +
  umap_scpas +
  umap_scpp +
  umap_scab +
  umap_degas +
  umap_lp_sgl +
  umap_pipet +
  umap_celltype +
  umap_subtype +
  plot_layout(ncol = 3)

ggplot2::ggsave(
  umaps,
  filename = "GSE9891_UMAP_screened.png",
  width = 15,
  height = 14,
  dpi = 400
)
