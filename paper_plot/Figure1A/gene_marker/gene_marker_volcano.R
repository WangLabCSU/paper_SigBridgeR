setwd(file.path(
  usethis::proj_path(),
  "paper_plot/Figure1A/gene_marker"
))

marker <- data.table::fread(
  "../../../4_positive_ctrl/luad/binary_deg_TCGA_LUAD.csv"
)

fc_threshold <- 1
fdr_threshold <- 0.05

plot_data <- marker |>
  dplyr::mutate(
    neg_log10_fdr = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
    direction = dplyr::case_when(
      adj.P.Val < fdr_threshold & logFC >= fc_threshold ~ "Upregulated",
      adj.P.Val < fdr_threshold & logFC <= -fc_threshold ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    direction = factor(
      direction,
      levels = c("Not significant", "Upregulated", "Downregulated")
    )
  )

label_data <- plot_data |>
  dplyr::filter(direction != "Not significant") |>
  dplyr::arrange(adj.P.Val, dplyr::desc(abs(logFC))) |>
  dplyr::slice_head(n = 20)

p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = logFC, y = neg_log10_fdr)) +
  ggplot2::geom_point(
    ggplot2::aes(color = direction),
    alpha = 0.7,
    size = 1.5
  ) +
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
  ggrepel::geom_text_repel(
    data = label_data,
    ggplot2::aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    seed = 5187
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Not significant" = "#CECECE",
      "Upregulated" = "#a02020",
      "Downregulated" = "#386c9b"
    )
  ) +
  ggplot2::labs(
    x = expression(log[2] * " fold change"),
    y = expression(-log[10] * " adjusted p-value"),
    color = NULL
  ) +
  cowplot::theme_cowplot() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 16),
    axis.title = ggplot2::element_text(size = 16),
    legend.title = ggplot2::element_text(size = 16),
    legend.text = ggplot2::element_text(size = 14)
  )

ggplot2::ggsave("gene_marker_volcano.png", p, width = 9, height = 7)
