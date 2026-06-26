setwd(file.path(usethis::proj_path(), "3_neg_ctrl_ucell/neg_viz"))

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)

ucell_score <- qs::qread(
  "../neg_ad_test/neg_score_nested_combined.qs",
  nthreads = 4L
) %>%
  purrr::imap(
    ~ {
      dtplyr::lazy_dt(.x) %>%
        tidyr::pivot_longer(
          cols = colnames(.x)[-1],
          names_to = "neg_sample",
          values_to = "ucell_mean_score"
        ) %>%
        dplyr::mutate(
          screen_method = gsub("_[^_]*$", "", cluster),
          data_name = .y,
          pheno = gsub("_.*", "", .y),
          bulk = if (stringr::str_detect(.y, "tcga")) {
            stringr::str_remove(
              .y,
              ".*_(?=[^_]+_[^_]+$)"
            )
          } else {
            stringr::str_remove(.y, ".*_(?=[^_]+$)")
          },
          sc = stringr::str_remove(.y, unique(bulk)) %>%
            stringr::str_remove(unique(pheno)) %>%
            stringr::str_remove_all("_"),
        ) %>%
        data.table::as.data.table()
    }
  ) %>%
  data.table::rbindlist()

ucell_score <- ucell_score[cluster != "SCIPAC_NA", ]

# 1. 确保 score_group 是因子，顺序固定
ucell_score$cluster <- factor(
  ucell_score$cluster,
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
    "LP_SGL_Neutral",
    "SCIPAC_Positive",
    "SCIPAC_Negative",
    "SCIPAC_Neutral"
  )
)
ucell_score$sc <- as.factor(
  ucell_score$sc
)
ucell_score$bulk <- as.factor(
  ucell_score$bulk
)
ucell_score$pheno <- as.factor(
  ucell_score$pheno
)

# 3. 获取实际的 x_group 水平（根据数据动态生成）
x_group_levels <- levels(ucell_score$cluster)
#   n_groups <- length(x_group_levels)

# 根据 x_group 名称自动分配颜色
palette <- sapply(x_group_levels, function(x) {
  if (grepl("Positive", x)) {
    return("#c24b4b")
  } else if (grepl("Negative", x)) {
    return("#5189bb")
  } else if (grepl("Neutral", x)) {
    return("#b6b6b6")
  } else if (grepl("Other", x)) {
    return("#b6b6b6")
  } else {
    return("#999999") # 默认颜色
  }
})
names(palette) <- x_group_levels


p <- ggplot2::ggplot(
  ucell_score,
  ggplot2::aes(
    x = cluster,
    y = ucell_mean_score,
    fill = cluster,
    color = cluster
  )
) +
  gghalves::geom_half_boxplot(
    side = "l",
    outlier.alpha = 0.2,
    outlier.size = 0.5,
    outlier.colour = "#cececeff",
    width = 0.65,
    alpha = 0.5,
    errorbar.length = 0.4,
    show.legend = FALSE
  ) +
  gghalves::geom_half_violin(
    side = "r",
    trim = FALSE,
    alpha = 0.5,
    show.legend = FALSE,
    width = 0.65,
    scale = "count"
  ) +
  #   ggplot2::geom_point(
  #     size = 0.3,
  #     alpha = 0.7,
  #     shape = 16,
  #     position = ggplot2::position_jitterdodge(
  #       jitter.width = 0.12,
  #       jitter.height = 0,
  #       dodge.width = 0.75
  #     )
  #   ) +
  # geom_signif(
  #   comparisons = list(c("A", "B"), c("A", "D")), # 设置需要比较的组
  #   map_signif_level = T, #是否使用星号显示
  #   test = t.test, ##计算方法
  #   y_position = c(25, 28), #图中横线位置设置
  #   tip_length = c(c(0.7, 0.3), c(0.8, 0.3)), #横线下方的竖线设置
  #   size = 1,
  #   color = "black"
  # )  +
  ggplot2::scale_fill_manual(values = palette, guide = "none") +
  ggplot2::scale_color_manual(values = palette, guide = "none") +
  cowplot::theme_cowplot() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 10
    ),
    axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
    axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      color = "#EEEEEE",
      linewidth = 0.3
    ),
    strip.text = ggplot2::element_text(face = "bold", size = 11),
    strip.background = ggplot2::element_rect(
      fill = "#EEEEEE",
      color = "white"
    )
  ) +
  ggplot2::labs(
    title = "Neg ctrl - UCell\n 100reps",
    x = "Screen Group",
    y = "UCell Mean Score"
  ) +
  ggplot2::facet_grid(
    sc + bulk ~ pheno,
    scales = "free"
  )


ggplot2::ggsave(
  filename = "ucell_100_test.png",
  plot = p,
  dpi = 400,
  width = 16,
  height = 24
)


p_hairtail <- ggplot2::ggplot(
  ucell_score,
  ggplot2::aes(
    x = cluster,
    y = ucell_mean_score,
    fill = cluster,
    color = cluster
  )
) +
  gghalves::geom_half_boxplot(
    side = "l",
    outlier.alpha = 0.2,
    outlier.size = 0.5,
    outlier.colour = "#cececeff",
    width = 0.65,
    alpha = 0.5,
    errorbar.length = 0.4,
    show.legend = FALSE
  ) +
  gghalves::geom_half_violin(
    side = "r",
    trim = FALSE,
    alpha = 0.5,
    show.legend = FALSE,
    width = 0.65,
    scale = "count"
  ) +
  #   ggplot2::geom_point(
  #     size = 0.3,
  #     alpha = 0.7,
  #     shape = 16,
  #     position = ggplot2::position_jitterdodge(
  #       jitter.width = 0.12,
  #       jitter.height = 0,
  #       dodge.width = 0.75
  #     )
  #   ) +
  # geom_signif(
  #   comparisons = list(c("A", "B"), c("A", "D")), # 设置需要比较的组
  #   map_signif_level = T, #是否使用星号显示
  #   test = t.test, ##计算方法
  #   y_position = c(25, 28), #图中横线位置设置
  #   tip_length = c(c(0.7, 0.3), c(0.8, 0.3)), #横线下方的竖线设置
  #   size = 1,
  #   color = "black"
  # )  +
  ggplot2::scale_fill_manual(values = palette, guide = "none") +
  ggplot2::scale_color_manual(values = palette, guide = "none") +
  cowplot::theme_cowplot() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 10
    ),
    axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
    axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      color = "#EEEEEE",
      linewidth = 0.3
    ),
    strip.text = ggplot2::element_text(face = "bold", size = 11),
    strip.background = ggplot2::element_rect(
      fill = "#EEEEEE",
      color = "white"
    )
  ) +
  ggplot2::labs(
    title = "Neg ctrl - UCell\n 100reps",
    x = "Screen Group",
    y = "UCell Mean Score"
  ) +
  ggplot2::facet_grid(
    bulk ~ sc + pheno,
    scales = "free"
  )


ggplot2::ggsave(
  filename = "ucell_100_test_hairtail.png",
  plot = p_hairtail,
  dpi = 400,
  width = 30,
  height = 16
)
