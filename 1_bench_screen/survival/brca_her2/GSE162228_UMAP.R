library(zeallot)
library(dplyr)

setwd(file.path(usethis::proj_path(), "1_bench_screen/survival/brca_her2"))

source("../../draw_umap.R")

data_path <- "/home/data/sigbridger/benchmark_data/brca/HER2"
bulk_name <- "GSE162228"
save_path <- file.path("plot", bulk_name)

dir.create(
  save_path,
  recursive = TRUE,
  showWarnings = FALSE
)

seurat_merged <- qs::qread(
  file.path(
    data_path,
    paste0("survival_her2_", bulk_name, "_merged_seurat.qs")
  ),
  nthreads = 8L
)

# seurat_her2_tumor <- readRDS(
#   "/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_HER2Tum.rds"
# )
# seurat_merged$is_tumor = ifelse(
#   colnames(seurat_merged) %in% rownames(seurat_her2_tumor@meta.data),
#   "TRUE",
#   "FALSE"
# )

# umap_cluster <- draw_umap(
#   seurat = seurat_merged,
#   group.by = "seurat_clusters",
#   title = "GSE161529 her2 seurat_clusters",
#   save_path = file.path(save_path, "GSE161529_her2_seurat_clusters_UMAP.png")
# )

# umap_tumor <- draw_umap(
#   seurat = seurat_merged,
#   group.by = "is_tumor",
#   cols = c("FALSE" = "#386c9b", "TRUE" = "#a02020"),
#   title = "GSE161529 her2 is tumor cell",
#   save_path = file.path(save_path, "GSE161529_her2_tumor_UMAP.png")
# )

c(
  umap_scissor,
  umap_scpas,
  umap_scipac,
  umap_scpp,
  umap_scab,
  umap_degas,
  umap_lp_sgl
  #   ,  umap_pipet
) %<-%
  purrr::map(
    c(
      "scissor",
      "scPAS",
      "SCIPAC",
      "scPP",
      'scAB',
      "DEGAS",
      "LP_SGL"
      #   ,      "PIPET"
    ),
    ~ draw_umap(
      seurat = seurat_merged,
      group.by = .x,
      label = FALSE,
      cols = c(
        "Other" = "#CECECE",
        "Neutral" = "#CECECE",
        "Positive" = "#c24b4b",
        "Negative" = "#5189bb"
      ),
      title = paste0("sc: GSE161529 her2\nbulk: ", bulk_name, "\nmethod: ", .x),
      save_path = file.path(
        save_path,
        paste0("GSE161529_her2_", bulk_name, "_", .x, "_UMAP.png")
      )
    ),
    .progress = "Drawing"
  )
