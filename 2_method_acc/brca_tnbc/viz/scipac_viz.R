setwd(file.path(usethis::proj_path(), "/2_method_acc/brca_tnbc"))

stats_dir <- "stats/"
method <- "scipac"

label_mats <- list.files(
  path = stats_dir,
  pattern = glue::glue("{method}.*\\.csv"),
  full.names = TRUE,
  ignore.case = TRUE
)
names(label_mats) <- basename(label_mats) %>%
  tools::file_path_sans_ext() %>%
  sub(".*_(?=[^_]*_[^_]*$)", "", ., perl = TRUE)

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
  if (name == "mat1_part1") {
    arg_samples <- data.frame(
      n_pc = sample(seq(10L, 100L, 10L), 50, replace = TRUE),
      resolution = sample(seq(0.1, 4L, 0.1), 50, replace = TRUE),
      ela_net_alpha = sample(seq(0.1, 1L, 0.1), 50, replace = TRUE)
    ) %>%
      dplyr::add_row(n_pc = 60L, resolution = 2L, ela_net_alpha = 0.4) # default
  } else if (name == "mat1_part2") {
    arg_samples <- data.frame(
      hvg = sample(seq(500L, 5000L, 500L), 50, replace = TRUE),
      bt_size = sample(seq(10L, 100L, 10L), 50, replace = TRUE),
      nfold = sample(seq(2L, 30L, 2L), 50, replace = TRUE)
    ) %>%
      dplyr::add_row(hvg = 1000L, bt_size = 50L, nfold = 10L) # default
  } else if (name == "mat2_part1") {
    arg_samples <- data.frame(
      n_pc = sample(seq(10L, 100L, 10L), 50, replace = TRUE),
      resolution = sample(seq(0.1, 4L, 0.1), 50, replace = TRUE),
      ela_net_alpha = sample(seq(0.1, 1L, 0.1), 50, replace = TRUE)
    ) %>%
      dplyr::add_row(n_pc = 60L, resolution = 2L, ela_net_alpha = 0.4) # default
  } else if (name == "mat2_part2") {
    arg_samples <- data.frame(
      hvg = sample(seq(500L, 5000L, 500L), 50, replace = TRUE),
      bt_size = sample(seq(10L, 100L, 10L), 50, replace = TRUE),
      nfold = sample(seq(2L, 30L, 2L), 50, replace = TRUE)
    ) %>%
      dplyr::add_row(hvg = 1000L, bt_size = 50L, nfold = 10L) # default
  }

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
  if (any(grepl("n_pc", colnames(data)))) {
    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = resolution, y = ela_net_alpha, fill = F1)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(size = n_pc, color = n_pc),
        alpha = 0.9,
        shape = 21,
        color = "black"
      ) +
      ggplot2::scale_fill_gradient(
        low = "white",
        high = "red",
        name = "F1"
      ) +
      # 1. 保证横轴按真实数值排序，并强制每 500 一个刻度
      ggplot2::scale_x_continuous(
        breaks = seq(0, 4, by = 0.5),
        labels = seq(0, 4, by = 0.5)
      ) +
      ggplot2::scale_y_continuous(
        breaks = seq(0, 1, by = 0.1),
        labels = seq(0, 1, by = 0.1)
      ) +
      ggplot2::labs(
        title = "Validation of the Screening Efficiency of SCIPAC under Random Parameters",
        subtitle = "x = resolution, y = ela_net_alpha",
        x = "resolution",
        y = "ela_net_alpha"
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
  } else {
    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = hvg, y = bt_size, fill = F1)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(size = nfold, color = nfold),
        alpha = 0.9,
        shape = 21,
        color = "black"
      ) +
      ggplot2::scale_fill_gradient(
        low = "white",
        high = "red",
        name = "F1"
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq(500, 5000, by = 500),
        labels = seq(500, 5000, by = 500)
      ) +
      ggplot2::scale_y_continuous(
        breaks = seq(10, 100, by = 10),
        labels = seq(10, 100, by = 10)
      ) +
      ggplot2::labs(
        title = "Validation of the Screening Efficiency of SCIPAC under Random Parameters",
        subtitle = "x = hvg, y = bt_size",
        x = "hvg",
        y = "bt_size"
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
  }

  ggplot2::ggsave(
    filename = save_path,
    plot = p,
    width = width,
    height = height,
    dpi = 400
  )

  p
}

mat1_p1 <- f1_bubble_heatmap(
  arg_samples_with_metrics$mat1_part1,
  save_path = glue::glue("viz/plot/{method}_acc1_part1.png")
)
mat1_p2 <- f1_bubble_heatmap(
  arg_samples_with_metrics$mat1_part2,
  save_path = glue::glue("viz/plot/{method}_acc1_part2.png")
)
mat2_p1 <- f1_bubble_heatmap(
  arg_samples_with_metrics$mat2_part1,
  save_path = glue::glue("viz/plot/{method}_acc2_part1.png")
)
mat2_p2 <- f1_bubble_heatmap(
  arg_samples_with_metrics$mat2_part2,
  save_path = glue::glue("viz/plot/{method}_acc2_part2.png")
)
