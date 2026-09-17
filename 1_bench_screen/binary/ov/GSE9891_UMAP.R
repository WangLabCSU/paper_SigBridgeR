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
if (
  !"SCIPAC" %in% colnames(seurat_merged[[]]) &&
    "sig" %in% colnames(seurat_merged[[]])
) {
  message("SCIPAC not in meta, use sig instead")
  # bug compatible
  seurat_merged$SCIPAC <- seurat_merged$sig
}

umap_cluster <- draw_umap(
  seurat = seurat_merged,
  group.by = "seurat_clusters",
  title = "Seurat Clusters",
  label.size = 4,
  save_path = file.path(save_path, "GSE165897_seurat_clusters_UMAP.png")
)

umap_tumor <- draw_umap(
  seurat = seurat_merged,
  group.by = "cell_type",
  label = FALSE,
  title = "Cell Type",
  save_path = file.path(save_path, "GSE165897_cell_type_UMAP.png")
)


method_cols <- c(
  "scissor",
  "scPAS",
  "SCIPAC",
  "scPP",
  "scAB",
  "DEGAS",
  "LP_SGL",
  "PIPET"
)
method_cols_missing <- setdiff(method_cols, colnames(seurat_merged[[]]))
if (length(method_cols_missing) > 0L) {
  cli::cli_alert_warning("Skip missing methods: {.val {method_cols_missing}}")
}
method_cols <- intersect(method_cols, colnames(seurat_merged[[]]))

umap_methods <-
  purrr::map(
    method_cols,
    ~ draw_umap(
      seurat = seurat_merged,
      group.by = .x,
      label = FALSE,
      cols = c(
        "Other" = "#CECECE",
        "Neutral" = "#CECECE",
        "Positive" = "#a02020",
        "Negative" = "#386c9b"
      ),
      title = glue::glue("Method: {.x}"),
      save_path = file.path(
        save_path,
        glue::glue("GSE165897_{bulk_name}_{.x}_UMAP.png")
      )
    ),
    .progress = "Drawing"
  )
