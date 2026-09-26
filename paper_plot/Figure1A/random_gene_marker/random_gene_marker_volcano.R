setwd(file.path(
  usethis::proj_path(),
  "paper_plot/Figure1A/random_gene_marker"
))

marker <- data.table::fread(
  "../../../4_positive_ctrl/luad/binary_deg_TCGA_LUAD.csv"
)
random_gene <- data.table::fread(
  "../../../3_negative_ctrl/luad/luad_random20_markers_100rep.csv"
)

fc_threshold <- 1
fdr_threshold <- 0.05
sample_columns <- paste0("Sample_", 1:20)

plot_data <- marker |>
  dplyr::mutate(
    neg_log10_fdr = -log10(pmax(adj.P.Val, .Machine$double.xmin))
  )

random_gene_data <- random_gene |>
  dplyr::select(dplyr::all_of(sample_columns)) |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "sample",
    values_to = "gene"
  ) |>
  dplyr::left_join(plot_data, by = "gene") |>
  dplyr::filter(!is.na(logFC)) |>
  dplyr::mutate(sample = factor(sample, levels = sample_columns))

sample_colors <- stats::setNames(
  grDevices::hcl.colors(length(sample_columns), palette = "Dark 3"),
  sample_columns
)

p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = logFC, y = neg_log10_fdr)) +
  ggplot2::geom_point(color = "#CECECE", alpha = 0.7, size = 1.5) +
  ggplot2::geom_vline(
    xintercept = c(-fc_threshold, fc_threshold),
    color = "#008B8B",
    linetype = "dashed"
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(fdr_threshold),
    color = "#008B8B",
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    data = random_gene_data,
    ggplot2::aes(color = sample),
    size = 1.8
  ) +
  ggplot2::scale_color_manual(values = sample_colors) +
  ggplot2::labs(
    x = expression(log[2] * " fold change"),
    y = expression(-log[10] * " adjusted p-value"),
    color = "Random sample"
  ) +
  cowplot::theme_cowplot() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 16),
    axis.title = ggplot2::element_text(size = 16),
    legend.title = ggplot2::element_text(size = 16),
    legend.text = ggplot2::element_text(size = 14)
  )

ggplot2::ggsave("random_gene_marker_volcano.png", p, width = 9, height = 7)
