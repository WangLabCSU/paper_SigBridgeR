DrawBulkPCA <- function(
  bulk,
  group,
  batch = NULL,
  show_plot = TRUE,
  ...
) {
  rlang::check_installed(c("ggplot2", "ggforce", "patchwork", "tibble"))

  pca <- stats::prcomp(t(bulk), scale. = TRUE)
  percent_var <- pca$sdev^2L / sum(pca$sdev^2L)

  pca_df <- tibble::tibble(
    sample = colnames(bulk),
    PC1 = pca$x[, 1L],
    PC2 = pca$x[, 2L],
    group = group,
    batch = batch
  )

  PC1 <- PC2 <- NULL
  p_pca <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = `PC1`, y = `PC2`, color = `group`)
  ) +
    ggplot2::geom_point(size = 3L, alpha = 0.8) +
    ggplot2::labs(
      x = paste0("PC1 (", round(percent_var[1L], 2L), "%)"),
      y = paste0("PC2 (", round(percent_var[2L], 2L), "%)")
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      legend.position = c(0.95, 0.95), # legend inside the top-right corner
      legend.justification = c(1L, 1L) # aligned to the top-right corner
    ) +
    ggforce::geom_mark_ellipse(
      ggplot2::aes(fill = group, group = group),
      alpha = 0.1,
      expand = ggplot2::unit(3L, "mm"),
      show.legend = FALSE
    ) +
    ggplot2::geom_hline(
      yintercept = 0L,
      color = "gray70",
      linetype = "dashed"
    ) +
    ggplot2::geom_vline(xintercept = 0L, color = "gray70", linetype = "dashed")

  if ("batch" %chin% colnames(pca_df)) {
    n_batch <- length(unique(pca_df$batch))

    chk::chk_lt(n_batch, 5L)

    p_pca <- p_pca +
      ggplot2::aes(shape = `batch`) +
      ggplot2::scale_shape_manual(
        values = c(
          16L,
          17L,
          18L,
          19L,
          20L
        )[seq_len(n_batch)]
      )
  }

  pc1_range <- range(pca_df$PC1)
  pc2_range <- range(pca_df$PC2)

  # create density plot - stroke color matches fill color
  density_x <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = PC1, fill = group, color = group)
  ) +
    ggplot2::geom_density(alpha = 0.7, bw = "nrd", adjust = 2L) +
    ggplot2::coord_cartesian(xlim = pc1_range) + # match the main plot x-axis range
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")

  density_y <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = PC2, fill = group, color = group)
  ) +
    ggplot2::geom_density(alpha = 0.7, trim = FALSE, bw = "nrd", adjust = 1L) +
    ggplot2::coord_cartesian(xlim = pc2_range) + # match the main plot y-axis range
    ggplot2::coord_flip() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")

  combined_ellipse <- density_x +
    patchwork::plot_spacer() +
    p_pca +
    density_y +
    patchwork::plot_layout(ncol = 2L, widths = c(4L, 1L), heights = c(1L, 4L)) +
    patchwork::plot_annotation(title = "Principal Component Analysis (PCA)")

  if (show_plot) {
    print(combined_ellipse)
  }

  invisible(combined_ellipse)
}
