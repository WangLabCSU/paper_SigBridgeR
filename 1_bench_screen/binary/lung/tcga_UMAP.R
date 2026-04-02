library(Seurat)
library(ggplot2)
library(patchwork)
library(zeallot)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


seurat_merged = qs::qread("tcga_luad_merged_seurat.qs", nthreads = 4)
set.seed(123)
my_color = randomcoloR::distinctColorPalette(
  k = length(unique(seurat_merged$seurat_clusters)),
  runTsne = TRUE
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
  umap_scab,
  # , umap_scpp,
  umap_degas,
  umap_lp_sgl,
  umap_pipet
) %<-%
  purrr::map(
    c(
      "scissor",
      "scPAS",
      "scAB",
      # , "scPP",
      "DEGAS",
      "LP_SGL",
      "PIPET"
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
  umap_degas +
  umap_lp_sgl +
  umap_pipet +
  umap_cnv +
  umap_tissue +
  plot_layout(ncol = 3)

ggplot2::ggsave(
  umaps,
  filename = "tcga_UMAP_screened.png",
  width = 16,
  height = 15,
  dpi = 400
)

# p = SigBridgeR::ScreenFractionPlot(
#   seurat_merged,
#   group_by = "Sample",
#   screen_type = c("scissor", "scPAS", "scAB")
# )
