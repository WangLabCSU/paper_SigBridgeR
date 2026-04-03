library(zeallot)
library(dplyr)

setwd(file.path(usethis::proj_path(), "1_bench_screen/binary/lung"))

source("../../draw_umap.R")

data_path <- "/home/data/sigbridger/benchmark_binary/lung"
bulk_name <- "TCGA_LUAD"
save_path <- file.path("plot", bulk_name)

dir.create(
  save_path,
  recursive = TRUE,
  showWarnings = FALSE
)

seurat_merged <- qs::qread(
  file.path(data_path, paste0("binary_lung_", bulk_name, "_merged_seurat.rqs")),
  nthreads = 8L
)

umap_cluster <- draw_umap(
  seurat = seurat_merged,
  group_by = "seurat_clusters",
  title = "GSE123902 seurat_clusters",
  save_path = file.path(save_path, "GSE123902_seurat_clusters_UMAP.png")
)

umap_tumor <- draw_umap(
  seurat = seurat_merged,
  group.by = "cnv_status",
  cols = c("normal" = "#386c9b", "tumor" = "#a02020"),
  save_path = file.path(save_path, "GSE123902_tumor_UMAP.png")
)


c(
  umap_scissor,
  umap_scpas,
  umap_scipac,
  umap_scpp,
  umap_scab,
  umap_degas,
  umap_lp_sgl,
  umap_pipet
) %<-%
  purrr::map(
    c(
      "scissor",
      "scPAS",
      "SCIPAC",
      "scPP",
      'scAB',
      "DEGAS",
      "LP_SGL",
      "PIPET"
    ),
    ~ draw_umap(
      seurat = seurat_merged,
      group.by = .x,
      cols = c(
        "Other" = "#CECECE",
        "Neutral" = "#CECECE",
        "Positive" = "#a02020",
        "Negative" = "#386c9b"
      ),
      title = paste0("sc: GSE123902\nbulk: ", bulk_name, "\nmethod: ", .x),
      save_path = file.path(
        save_path,
        paste0("GSE123902_", bulk_name, "_", .x, "_UMAP.png")
      )
    )
  )
