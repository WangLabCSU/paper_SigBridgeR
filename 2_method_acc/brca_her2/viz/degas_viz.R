setwd(file.path(usethis::proj_path(), "/2_method_acc/brca_her2"))

stats_dir <- "stats/"
method <- "degas"

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
  set.seed(123)
  part_name <- gsub("mat._", "", name)
  if (part_name == "part1") {
    arg_samples <- data.frame(
      arch = sample(c("DenseNet", "Standard"), 50, replace = TRUE),
      ff_depth = sample(2:10, 50, replace = TRUE),
      bag_depth = sample(3:10, 50, replace = TRUE)
    ) %>%
      dplyr::add_row(arch = "DenseNet", ff_depth = 3, bag_depth = 5)
  } else if (part_name == "part2") {
    arg_samples <- data.frame(
      lamb1 = sample(2:10, 50, replace = TRUE),
      lamb2 = sample(2:10, 50, replace = TRUE),
      lamb3 = sample(2:10, 50, replace = TRUE)
    ) %>%
      dplyr::add_row(lamb1 = 3, lamb2 = 3, lamb3 = 3)
  } else if (part_name == "part3") {
    arg_samples <- data.frame(
      scbatch_sz = sample(seq(50, 500, 10), 50, replace = TRUE),
      patbatch_sz = sample(seq(25, 100, 5), 50, replace = TRUE),
      hidden_feats = sample(seq(25, 100, 5), 50, replace = TRUE),
      do_prc = sample(seq(0.1, 0.9, 0.1), 50, replace = TRUE)
    ) %>%
      dplyr::add_row(
        scbatch_sz = 200,
        patbatch_sz = 50,
        hidden_feats = 50,
        do_prc = 0.5
      )
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


# -----------------------------------------------------------------------------------------------
# * viz

f1_bubble_heatmap <- function(
  data,
  save_path = "viz/plot",
  width = 10,
  height = 8,
  ...
) {
  if (any(grepl("arch", colnames(data)))) {
    # arg_samples1
    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = ff_depth, y = bag_depth, fill = F1)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(size = arch, color = arch),
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
        breaks = seq(2, 10, by = 1),
        labels = seq(2, 10, by = 1)
      ) +
      ggplot2::scale_y_continuous(
        breaks = seq(3, 10, by = 1),
        labels = seq(3, 10, by = 1)
      ) +
      ggplot2::labs(
        title = "Validation of the Screening Efficiency of DEGAS under Random Parameters",
        subtitle = "x = ff_depth, y = bag_depth",
        x = "ff_depth",
        y = "bag_depth"
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
  } else if (any(grepl("lamb1", colnames(data)))) {
    # arg_samples2
    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = lamb1, y = lamb2, fill = F1)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(size = lamb3, color = lamb3),
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
        breaks = seq(2, 10, by = 1),
        labels = seq(2, 10, by = 1)
      ) +
      ggplot2::scale_y_continuous(
        breaks = seq(3, 10, by = 1),
        labels = seq(3, 10, by = 1)
      ) +
      ggplot2::labs(
        title = "Validation of the Screening Efficiency of DEGAS under Random Parameters",
        subtitle = "x = lamb1, y = lamb2",
        x = "lamb1",
        y = "lamb2"
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
  if (any(grepl("do_prc", colnames(data)))) {
    # arg_samples3
    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = scbatch_sz, y = patbatch_sz, fill = F1)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(size = hidden_feats, color = hidden_feats),
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
        breaks = seq(50, 500, by = 50),
        labels = seq(50, 500, by = 50)
      ) +
      ggplot2::scale_y_continuous(
        breaks = seq(25, 100, by = 25),
        labels = seq(25, 100, by = 25)
      ) +
      ggplot2::labs(
        title = "Validation of the Screening Efficiency of SCIPAC under Random Parameters",
        subtitle = "x = scbatch_sz, y = patbatch_sz",
        x = "scbatch_sz",
        y = "patbatch_sz"
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
      ) +
      ggplot2::facet_wrap(~do_prc)
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

purrr::iwalk(
  .x = arg_samples_with_metrics,
  .f = \(data, name) {
    mat_index <- gsub("mat(.).*", "\\1", name)
    part_name <- gsub("mat._", "", name)

    f1_bubble_heatmap(
      data = data,
      save_path = glue::glue("viz/plot/{method}_acc{mat_index}_{part_name}.png")
    )
  },
  .progress = "Drawing"
)
