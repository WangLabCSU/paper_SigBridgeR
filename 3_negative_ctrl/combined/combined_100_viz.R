setwd(file.path(usethis::proj_path(), "3_negative_ctrl/combined"))

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)

data_paths <- list.files(
  path = "..",
  pattern = ".stat\\.csv",
  recursive = TRUE,
  full.names = TRUE
)
bulk <- dirname(data_paths) %>% basename()
sc_type <- dirname(dirname(data_paths)) %>% basename() %>% toupper()
pheno_type <- gsub(".*(survival|binary).*", "\\1", data_paths)

names(data_paths) <- paste0(pheno_type, "_", sc_type, "_", bulk)

data_loaded <- purrr::imap(
  .x = data_paths,
  .f = function(path, name) {
    dt <- data.table::fread(path)
    dt[, data_name := name] %>%
      tidyr::separate(
        col = "data_name",
        remove = FALSE,
        into = c("pheno_type", "sc_type", "bulk")
      )
  }
)

data_combined <- dplyr::bind_rows(data_loaded)

# 1. 确保 score_group 是因子，顺序固定
data_combined$score_group <- as.factor(
  data_combined$score_group
)

# 2. 创建组合分组标签 (用于 x 轴)
data_combined <- dplyr::mutate(
  data_combined,
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
        "LP_SGL_Neutral",
        "SCIPAC_Positive",
        "SCIPAC_Negative",
        "SCIPAC_Neutral"
      )
    )
  )

# 3. 获取实际的 x_group 水平（根据数据动态生成）
x_group_levels <- levels(data_combined$x_group)
#   n_groups <- length(x_group_levels)

# 4. 创建调色板（如果未提供）

# 根据 x_group 名称自动分配颜色
palette <- sapply(x_group_levels, function(x) {
  if (grepl("Positive", x)) {
    return("#d65456ff")
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


p <- ggplot2::ggplot(
  data_combined,
  ggplot2::aes(x = x_group, y = mean_score, fill = x_group, color = x_group)
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
    title = "Neg ctrl",
    x = "Screen Group",
    y = "ssGSEA Mean Score"
  ) +
  ggplot2::facet_grid(
    sc_type + bulk ~ pheno_type,
    scales = "free"
  )


ggplot2::ggsave(
  filename = "combined_100_test.png",
  plot = p,
  dpi = 400,
  width = 16,
  height = 24
)
