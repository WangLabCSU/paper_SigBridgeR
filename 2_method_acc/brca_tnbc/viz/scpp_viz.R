setwd(file.path(usethis::proj_path(), "/2_method_acc/brca_tnbc"))

stats_dir <- "stats/"
method <- "scpp"

label_mats <- list.files(
  path = stats_dir,
  pattern = glue::glue("{method}.*\\.csv"),
  full.names = TRUE,
  ignore.case = TRUE
)
names(label_mats) <- basename(label_mats) %>%
  tools::file_path_sans_ext() %>%
  gsub(".*_", "", .)

label_mats_loaded <- lapply(label_mats, data.table::fread)


# * benchmark label
data_dir <- "/home/data/sigbridger/benchmark_data/brca"
sc_data <- qs::qread(file.path(data_dir, "seurat_tnbc.qs"), nthreads = 8L)
seurat_tumor <- readRDS(
  "/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_TNBCTum.rds"
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
    prob = sample(seq(0, 0.5, 0.01), 50, replace = TRUE),
    Log2FC_cutoff = round(runif(50), 3)
  ) %>%
    dplyr::add_row(prob = 0.2, Log2FC_cutoff = 0.585) # default

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


# --------------------------------------------------------------------------------------------------------------------------
# * viz

f1_bubble_heatmap <- function(
  data,
  save_path = "viz/plot",
  width = 10,
  height = 8,
  ...
) {
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = prob, y = Log2FC_cutoff, fill = F1)
  ) +
    ggplot2::geom_point(size = 6, alpha = 0.9, shape = 21, color = "black") +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "red",
      name = "F1"
    ) +
    # 1. 保证横轴按真实数值排序，并强制每 500 一个刻度
    ggplot2::scale_x_continuous(
      breaks = seq(0, 1, by = 0.05),
      labels = seq(0, 1, by = 0.05)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, by = 0.1),
      labels = seq(0, 1, by = 0.1)
    ) +
    ggplot2::labs(
      title = "Validation of the Screening Efficiency of scPP under Random Parameters",
      subtitle = "x = prob, y = log2FC_cutoff",
      x = "prob",
      y = "log2FC_cutoff"
    ) +
    ggplot2::theme_minimal(base_size = 14) + # 全局字体基准
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 12), # 2. 轴文字放大
      axis.title = ggplot2::element_text(size = 13),
      # 3. x 轴 45° 倾斜
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1),
      # 4. 图例文字放大
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 13)
    )

  ggplot2::ggsave(
    filename = save_path,
    plot = p,
    width = width,
    height = height,
    dpi = 400
  )

  p
}

mat1_p <- f1_bubble_heatmap(
  arg_samples_with_metrics$mat1,
  save_path = "viz/plot/scpp_acc1.png"
)
mat2_p <- f1_bubble_heatmap(
  arg_samples_with_metrics$mat2,
  save_path = "viz/plot/scpp_acc2.png"
)
