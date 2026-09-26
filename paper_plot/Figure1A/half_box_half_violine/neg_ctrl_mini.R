setwd(file.path(
  usethis::proj_path(),
  "paper_plot/Figure1A/half_box_half_violine"
))

library(dplyr)
library(ggplot2)
library(gghalves)

data <- qs::qread(
  "../../../3_negative_ctrl/luad/luad_Sample_100_ssgsea_score.qs"
)
seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/TCGA-LUAD/luad_tcga_scissor.qs",
  nthreads = 8L
)
label <- seurat$scRNA_data$scissor

stopifnot(nrow(data) == length(label))

plot_data <- tibble(
  label = factor(label, levels = c("Positive", "Negative", "Neutral")),
  mean_score = rowMeans(data, na.rm = TRUE)
) |>
  filter(!is.na(label), is.finite(mean_score))

palette <- c(
  Positive = "#a02020",
  Negative = "#386c9b",
  Neutral = "#CECECE"
)

p <- ggplot(
  plot_data,
  aes(x = label, y = mean_score, fill = label, color = label)
) +
  geom_half_boxplot(
    side = "l",
    outlier.alpha = 0.2,
    outlier.size = 0.5,
    outlier.colour = "#CECECE",
    width = 0.65,
    alpha = 0.5,
    errorbar.length = 0.4,
    show.legend = FALSE
  ) +
  geom_half_violin(
    side = "r",
    trim = FALSE,
    alpha = 0.5,
    width = 0.65,
    scale = "area",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = palette, guide = "none") +
  scale_color_manual(values = palette, guide = "none") +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = 0.3)
  ) +
  labs(
    x = "Scissor classification",
    y = "Mean negative-control ssGSEA score"
  )

ggsave(
  filename = "neg_ctrl_mini.png",
  plot = p,
  dpi = 400,
  width = 5,
  height = 4
)
