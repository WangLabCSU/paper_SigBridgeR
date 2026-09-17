library(zeallot)
library(dplyr)

setwd(file.path(usethis::proj_path(), "1_bench_screen/binary/fshd"))

source("../../draw_umap.R")

data_path <- "/home/data/sigbridger/benchmark_binary/fshd"
bulk_name <- "GSE56787"
save_path <- file.path("plot", bulk_name)

dir.create(
  save_path,
  recursive = TRUE,
  showWarnings = FALSE
)

seurat_merged <- qs::qread(
  file.path(data_path, paste0("binary_fshd_", bulk_name, "_merged_seurat.qs")),
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
  save_path = file.path(
    save_path,
    glue::glue("GSE122873_seurat_clusters_UMAP.png")
  )
)

# CTRL / FSHD1 / FSHD2
umap_condition <- draw_umap(
  seurat = seurat_merged,
  group.by = "condition",
  label.size = 4,
  cols = c("CTRL" = "#386c9b", "FSHD1" = "#a02020", "FSHD2" = "#e0822e"),
  title = "Condition",
  save_path = file.path(
    save_path,
    glue::glue("GSE122873_condition_UMAP.png")
  )
)

# healthy control vs FSHD
umap_disease_state <- draw_umap(
  seurat = seurat_merged,
  group.by = "ch1_disease state",
  label = FALSE,
  cols = c(
    "healthy control" = "#386c9b",
    "Facioscapulohumeral muscular dystrophy" = "#a02020"
  ),
  title = "Disease State",
  save_path = file.path(
    save_path,
    glue::glue("GSE122873_disease_state_UMAP.png")
  ),
  width = 12
)

c(
  umap_sample,
  umap_sex
) %<-%
  purrr::map(
    c("sample", "ch1_Sex"),
    ~ draw_umap(
      seurat = seurat_merged,
      group.by = .x,
      label.size = 4,
      title = glue::glue("{.x}"),
      save_path = file.path(
        save_path,
        glue::glue("GSE122873_{.x}_UMAP.png")
      )
    ),
    .progress = "Drawing"
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
        glue::glue("GSE122873_{bulk_name}_{.x}_UMAP.png")
      )
    ),
    .progress = "Drawing"
  )
