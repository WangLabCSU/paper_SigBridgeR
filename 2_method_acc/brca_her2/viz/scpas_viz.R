setwd(file.path(usethis::proj_path(), "/2_method_acc/brca_her2"))

stats_dir <- "stats/"
method <- "scpas"

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
  if (name == "mat1") {
    set.seed(123)
  } else if (name == "mat2") {
    set.seed(42568)
  } else {
    cli::cli_abort("Unknown name: {.val {name}}")
  }
  arg_samples <- data.frame(
    "nfeature" = sample(seq(500, 5000, by = 100), 50, replace = TRUE), # 第1维
    "imputation" = sample(c("None", "KNN", "ALRA"), 50, replace = TRUE), # 第2维
    "independent" = sample(c(FALSE, TRUE), 50, replace = TRUE) # 第3维
  ) %>%
    dplyr::add_row(nfeature = 2000, imputation = "None", independent = TRUE)

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
  data <- tidyr::unite(
    data = data,
    col = "y_axis",
    imputation,
    independent,
    sep = " | "
  )

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = nfeature, y = y_axis, fill = F1)
  ) +
    ggplot2::geom_point(size = 6, alpha = 0.9, shape = 21, color = "black") +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "red",
      name = "F1"
    ) +
    # 1. 保证横轴按真实数值排序，并强制每 500 一个刻度
    ggplot2::scale_x_continuous(
      breaks = seq(500, 5000, by = 250),
      labels = seq(500, 5000, by = 250)
    ) +
    ggplot2::labs(
      title = "Validation of the Screening Efficiency of scPAS under Random Parameters",
      subtitle = "x = nfeature, y = imputation | independent",
      x = "Number of features",
      y = NULL
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
  save_path = "viz/plot/scpas_acc1.png"
)
mat2_p <- f1_bubble_heatmap(
  arg_samples_with_metrics$mat2,
  save_path = "viz/plot/scpas_acc2.png"
)
