setwd(usethis::proj_path())

library(dplyr)

PlotGroupMeanScores <- function(
  data,
  palette = NULL,
  point_size = 2.5,
  outlier_size = 3.5,
  outlier_shape = 17,
  box_width = 0.6,
  dodge_width = 0.7,
  title = "Mean Score Distribution by Group",
  subtitle = "GSE42568",
  x_label = "Meta Column & Label",
  y_label = "Mean Score",
  angle = 45
) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggsignif)

  # 1. 确保 score_group 是因子，顺序固定
  score_levels <- unique(data$score_group)
  data$score_group <- factor(data$score_group, levels = score_levels)

  # 2. 创建组合分组标签 (用于 x 轴)
  data <- dplyr::mutate(
    data,
    x_group = interaction(meta_column, label_value, sep = "_", drop = TRUE)
  ) %>%
    dplyr::mutate(
      x_group = factor(
        x_group,
        levels = c(
          "scAB_Positive",
          "scAB_Other",
          "DEGAS_Positive",
          "DEGAS_Other",
          "PIPET_Positive",
          "PIPET_Negative",
          "PIPET_Neutral",
          "scPAS_Positive",
          "scPAS_Negative",
          "scPAS_Neutral",
          "scPP_Positive",
          "scPP_Negative",
          "scPP_Neutral",
          "scissor_Positive",
          "scissor_Negative",
          "scissor_Neutral",
          "LP_SGL_Positive",
          "LP_SGL_Negative",
          "LP_SGL_Neutral"
        )
      )
    )

  # 3. 获取实际的 x_group 水平（根据数据动态生成）
  x_group_levels <- levels(data$x_group)
  n_groups <- length(x_group_levels)

  # 4. 创建调色板（如果未提供）
  if (is.null(palette)) {
    # 根据 x_group 名称自动分配颜色
    palette <- sapply(x_group_levels, function(x) {
      if (grepl("Positive", x)) {
        return("#ff3333")
      } else if (grepl("Negative", x)) {
        return("#386c9b")
      } else if (grepl("Neutral", x)) {
        return("#b6b6b6")
      } else if (grepl("Other", x)) {
        return("#b6b6b6")
      } else {
        return("#999999") # 默认颜色
      }
    })
    names(palette) <- x_group_levels
  }

  p <- ggplot(
    data,
    aes(x = x_group, y = mean_score, fill = x_group, color = x_group)
  ) +
    geom_boxplot(
      width = box_width,
      outlier.size = outlier_size,
      outlier.shape = outlier_shape,
      outlier.color = "black",
      outlier.stroke = 0.5,
      outlier.alpha = 0.9,
      alpha = 0.7,
      position = position_dodge(width = dodge_width)
    ) +
    geom_point(
      size = point_size,
      alpha = 0.7,
      shape = 16,
      position = position_jitterdodge(
        jitter.width = 0.12,
        jitter.height = 0,
        dodge.width = dodge_width
      )
    ) +
    # geom_signif(
    #   comparisons = list(c("A", "B"), c("A", "D")), # 设置需要比较的组
    #   map_signif_level = T, #是否使用星号显示
    #   test = t.test, ##计算方法
    #   y_position = c(25, 28), #图中横线位置设置
    #   tip_length = c(c(0.7, 0.3), c(0.8, 0.3)), #横线下方的竖线设置
    #   size = 1,
    #   color = "black"
    # )  +
    scale_fill_manual(values = palette, guide = "none") +
    scale_color_manual(values = palette, guide = "none") +
    cowplot::theme_cowplot() +
    theme(
      axis.text.x = element_text(
        angle = angle,
        hjust = 1,
        vjust = 1,
        size = 10
      ),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5)
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label
    )

  return(p)
}

mean_scores <- data.table::fread(
  "Tmp/ssGSEA_negative_compare/survival/brca/her2/GSE42568/her2_sur_100reps_neg_ctrl_stat.csv"
)

mean_scores <- tidyr::unite(
  mean_scores,
  "ad_group",
  meta_column,
  label_value,
  sep = "_",
  remove = FALSE
)


p <- PlotGroupMeanScores(
  mean_scores,
  title = "Negative Control - Random ssGSEA Mean Score Distribution by Group",
  subtitle = "GSE42568 survival",
  x_label = "Screening Group",
  y_label = "Mean Score of ssGSEA",
)
p

ggplot2::ggsave(
  "Tmp/ssGSEA_negative_compare/survival/brca/her2/GSE42568/her2_sur_100reps_neg_ctrl_stat.png",
  p,
  width = 8,
  height = 6,
  dpi = 400
)
