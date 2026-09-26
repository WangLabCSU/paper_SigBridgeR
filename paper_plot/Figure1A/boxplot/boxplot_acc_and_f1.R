setwd(file.path(
  usethis::proj_path(),
  "paper_plot/Figure1A/boxplot"
))

data <- data.table::fread("../../../2_method_acc/viz/combined.csv")

plot_data <- data |>
  dplyr::select(F1, Accuracy) |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    metric = factor(metric, levels = c("F1", "Accuracy"))
  )

p <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = metric, y = value, fill = metric)
) +
  ggplot2::geom_boxplot(
    width = 0.55,
    color = "black",
    outlier.shape = 16,
    outlier.size = 1.5,
    fatten = 1.5
  ) +
  ggplot2::stat_boxplot(
    geom = "errorbar",
    width = 0.25,
    color = "black",
    linewidth = 0.6
  ) +
  ggplot2::stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.55,
    color = "black",
    linewidth = 0.8,
    fill = NA
  ) +
  ggplot2::scale_fill_manual(
    values = c("F1" = "#B8D9A0", "Accuracy" = "#9BBBD9")
  ) +
  ggplot2::labs(x = NULL, y = "Metric value") +
  cowplot::theme_cowplot() +
  ggplot2::theme(
    legend.position = "none",
    axis.title = ggplot2::element_text(size = 16),
    axis.text.x = ggplot2::element_text(size = 16),
    axis.text.y = ggplot2::element_text(size = 14)
  )

ggplot2::ggsave("boxplot_acc_and_f1.png", p, width = 5, height = 5, dpi = 400)
