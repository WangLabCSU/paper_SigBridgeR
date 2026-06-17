# ! HER2
# ! Run this after `compute_hallmarks_score.R`

library(dplyr)

# * edit this variables
tissue <- "her2"
survival_seurats_dir <- "/home/data/sigbridger/benchmark_data/brca/HER2"
binary_seurats_dir <- "/home/data/sigbridger/benchmark_binary/brca/HER2"

gsea_res <- list.files(
  "/home/data/sigbridger/GSEA",
  pattern = "gsea",
  full.names = TRUE
)
names(gsea_res) <- basename(gsea_res) %>%
  tools::file_path_sans_ext() %>%
  gsub("_.*", "", .)

gsea_score <- data.table::fread(gsea_res[[tissue]][[1]])

survival_seurats <- list.files(
  survival_seurats_dir,
  pattern = "survival",
  full.names = TRUE
)
names(survival_seurats) <- basename(survival_seurats) %>%
  tools::file_path_sans_ext()
binary_seurats <- list.files(
  binary_seurats_dir,
  pattern = "binary",
  full.names = TRUE
)
names(binary_seurats) <- basename(binary_seurats) %>%
  tools::file_path_sans_ext()

seurats <- c(survival_seurats, binary_seurats)

cli::cli_alert_info("Aggregating meta data with GSEA score")

methods <- c(
  "scissor",
  "scAB",
  "scPAS",
  "scPP",
  "DEGAS",
  "PIPET",
  "SCIPAC",
  "LP_SGL"
)
output_dir <- "/home/data/sigbridger/GSEA"

for (i in seq_along(seurats)) {
  seurat <- qs::qread(seurats[[i]], nthreads = 8L)
  existing_cols <- colnames(seurat[[]])
  filtered_cols <- existing_cols[existing_cols %in% methods]
  meta <- seurat[[]] %>% dplyr::select(dplyr::all_of(filtered_cols))
  seurat_name <- names(seurats)[[i]]

  cli::cli_alert_info("{seurat_name}: {.val {filtered_cols}}")

  gsea_score_combined <- cbind(gsea_score, meta)

  data.table::fwrite(
    gsea_score_combined,
    file.path(output_dir, paste0(seurat_name, "gsea_score.csv"))
  )
  gc(verbose = FALSE)
}
