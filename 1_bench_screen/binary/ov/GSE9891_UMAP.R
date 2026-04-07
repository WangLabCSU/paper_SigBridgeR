library(zeallot)
library(dplyr)

setwd(file.path(usethis::proj_path(), "1_bench_screen/binary/ov"))

source("../../draw_umap.R")

data_path <- "/home/data/sigbridger/benchmark_binary/ov"
bulk_name <- "GSE9891"
save_path <- file.path("plot", bulk_name)

dir.create(
  save_path,
  recursive = TRUE,
  showWarnings = FALSE
)

seurat_merged <- qs::qread(
  file.path(data_path, paste0("binary_ov_", bulk_name, "_merged_seurat.qs")),
  nthreads = 8L
)

umap_cluster <- draw_umap(
  seurat = seurat_merged,
  group.by = "seurat_clusters",
  title = "GSE165897 seurat_clusters",
  label.size = 4,
  save_path = file.path(save_path, "GSE165897_seurat_clusters_UMAP.png")
)

umap_tumor <- draw_umap(
  seurat = seurat_merged,
  group.by = "cell_type",
  label = FALSE,
  cols = c("Stromal" = "#70b641", "Immune" = "#b3bd30", "EOC" = "#a02020"),
  save_path = file.path(save_path, "GSE165897_tumor_UMAP.png")
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
      label = FALSE,
      cols = c(
        "Other" = "#CECECE",
        "Neutral" = "#CECECE",
        "Positive" = "#c24b4b",
        "Negative" = "#5189bb"
      ),
      title = paste0("sc: GSE165897\nbulk: ", bulk_name, "\nmethod: ", .x),
      save_path = file.path(
        save_path,
        paste0("GSE165897_", bulk_name, "_", .x, "_UMAP.png")
      )
    ),
    .progress = "Drawing"
  )
