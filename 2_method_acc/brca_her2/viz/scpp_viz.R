setwd(file.path(usethis::proj_path(), "/2_method_acc/brca_her2"))

stats_dir <- "stats/"
method <- "scpp"

label_mats <- list.files(
  path = stats_dir,
  pattern = method,
  full.names = TRUE,
  ignore.case = TRUE
)
names(label_mats) <- basename(label_mats) %>%
  tools::file_path_sans_ext() %>%
  gsub(".*_", "", .)

label_mats_loaded <- lapply(label_mats, data.table::fread)


# * benchmark label
data_dir <- "/home/data/sigbridger/benchmark_data/brca"
sc_data <- qs::qread(file.path(data_dir, "seurat_her2.qs"), nthreads = 8L)
seurat_tumor <- readRDS(
  "/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_HER2Tum.rds"
)
tumor_cells <- rownames(seurat_tumor@meta.data)
benchmark_label <- colnames(sc_data) %in% tumor_cells

# * add bench col to data
label_mats_loaded <- lapply(label_mats_loaded, function(x) {
  x[, benchmark := benchmark_label]
})

# * get the function
source("../compute_metrics.R")

# * compute metrics
metrics <- lapply(label_mats_loaded, function(dt) {
  compute_metrics(dt)
})

# * add arg_samples to metrics and save this result
arg_samples_with_metrics <- purrr::imap(metrics, \(dt, name) {
  set.seed(12345)

  arg_samples <- data.frame(
    prob = "NULL",
    Log2FC_cutoff = round(runif(50), 3)
  ) %>%
    dplyr::add_row(prob = "0.2", Log2FC_cutoff = 0.585) # default

  col_names <- rownames(dt)
  t_dt <- t(dt)
  colnames(t_dt) <- col_names
  cbind(arg_samples, t_dt)
})

purrr::iwalk(arg_samples_with_metrics, function(dt, name) {
  index <- gsub("mat", "", name)

  data.table::fwrite(
    dt,
    file.path(
      "arg_samples",
      glue::glue("{tolower(method)}_arg_samples{index}.csv")
    )
  )
  cli::cli_alert_success("{name} saved")
})
