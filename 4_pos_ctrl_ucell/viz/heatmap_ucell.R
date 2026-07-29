setwd(file.path(usethis::proj_path(), "4_pos_ctrl_ucell/viz/"))

library(data.table)

data_combined <- data.table::fread("../diff_test/all_diff_df.csv")

# ? all method comparisons, used to complete missing combinations
method_comparisons <- expand.grid(
  bulk = unique(data_combined$bulk),
  comparison = purrr::map_vec(
    c("scissor", "scAB", "scPAS", "scPP", "SCIPAC", "DEGAS", "LP_SGL", "PIPET"),
    ~ paste(.x, "Positive vs Rest", sep = "_")
  )
)

# ? generate triangle coordinates
triangle <- function(pairs, type = "up") {
  # 默认的上三角坐标基
  x = c(0, 0, 1)
  y = c(0, 1, 1)
  # 下三角的坐标基
  if (type == "lower") {
    x = c(0, 1, 1)
    y = c(0, 0, 1)
  }
  # 生成三角矩阵
  mat = do.call(
    rbind,
    apply(pairs, 1, function(row) {
      a = row[1]
      b = row[2]
      data.frame(
        x = x + a,
        y = y + b,
        group = paste(a, b, sep = "-")
      )
    })
  )
  return(mat)
}

# ? generate triangle data for ggplot2
triangle_data <- function(data, row = 1, col = 2) {
  # 保留原始因子水平顺序
  row_levels <- unique(data[[row]])
  col_levels <- unique(data[[col]])

  rows <- length(row_levels)
  cols <- length(col_levels)
  pairs <- merge(1:rows, 1:cols)

  # 获取上三角坐标
  upper <- triangle(pairs)
  colnames(upper) <- c(paste0("upper.", colnames(upper)[1:2]), "group")

  # 获取下三角坐标
  lower <- triangle(pairs, type = "lower")[1:2]
  colnames(lower) <- paste0("lower.", colnames(lower))

  # 合并坐标
  upper_lower <- bind_cols(upper, lower)
  pairs <- tidyr::unite(data = pairs, col = "group", x, y, sep = "-")

  # **关键修改**：创建映射，使用原始顺序而非转换后的数值
  data$row_numeric <- as.numeric(factor(data[[row]], levels = row_levels))
  data$col_numeric <- as.numeric(factor(data[[col]], levels = col_levels))
  data$group <- paste(data$row_numeric, data$col_numeric, sep = "-")

  # 根据分组信息将坐标连接到数据中
  data <- data %>%
    right_join(upper_lower, by = "group")

  return(data)
}

# ? complete missing combinations with NA
# 在调用 triangle_data 之前
plot_df <- data_combined %>%
  dtplyr::lazy_dt() %>%
  dplyr::mutate(
    UCell_type = gsub(".*(pos|neg).*", "\\1", score),
    comparison = paste0(screen_method, "_Positive vs Rest")
  ) %>%
  dplyr::mutate(
    # 将 bulk 和 comparison 转换为因子，并指定你想要的顺序
    comparison = factor(
      comparison,
      levels = unique(sort(comparison))
    ),
    bulk = factor(bulk, levels = unique(sort(bulk))),
    sc = factor(sc, levels = unique(sort(sc))),
    UCell_type = factor(UCell_type, levels = unique(sort(UCell_type)))
  ) %>%
  dplyr::arrange(bulk, comparison) %>%
  data.table::as.data.table()

plot_df2 <- triangle_data(plot_df, row = "bulk", col = "comparison")

plot_df2 <- plot_df2[!is.na(bulk)] %>%
  dtplyr::lazy_dt() %>%
  dplyr::mutate(
    neg_log10_p = -log10(p_value + 1e-300), # 300 是上限
    label = case_when(
      p_value < 1e-3 ~ "***",
      p_value < 1e-2 ~ "**",
      p_value < 0.05 ~ "*",
      p_value >= 0.05 ~ "NS",
      TRUE ~ "NA"
    )
  ) %>%
  data.table::as.data.table()


# ? significance label position
centers <- plot_df2 %>%
  dtplyr::lazy_dt() %>%
  dplyr::group_by(
    bulk,
    comparison,
    pheno,
    UCell_type,
    sc
  ) %>%
  dplyr::summarise(
    x_center = mean(lower.x),
    y_center = mean(lower.y),
    label = dplyr::first(label),
    .groups = "drop"
  ) %>%
  dplyr::filter(label != "" & !is.na(label)) %>%
  data.table::as.data.table()

plot_df2[sc == "her2", sc := "brca her2"]
centers[sc == "her2", sc := "brca her2"]
plot_df2[,
  bulk := toupper(bulk)
][,
  bulk := factor(bulk, levels = unique(sort(bulk)))
]
centers[,
  bulk := toupper(bulk)
][,
  bulk := factor(bulk, levels = unique(sort(bulk)))
]
plot_df2[,
  UCell_type := data.table::fcase(
    UCell_type == "pos" , "pos_ssGSEA" ,
    UCell_type == "neg" , "neg_ssGSEA"
  )
]
centers[,
  UCell_type := data.table::fcase(
    UCell_type == "pos" , "pos_ssGSEA" ,
    UCell_type == "neg" , "neg_ssGSEA"
  )
]

x_breaks <- seq_along(unique(plot_df2[["bulk"]])) + 0.5
x_labels <- sort(unique(plot_df2[["bulk"]]))

y_breaks <- seq_along(unique(plot_df2[["comparison"]])) + 0.5
y_labels <- sort(unique(plot_df2[["comparison"]]))


p <- ggplot2::ggplot(plot_df2) +
  ggplot2::geom_polygon(
    ggplot2::aes(upper.x, upper.y, fill = diff, group = group),
    colour = "grey",
    linewidth = 0.1
  ) +
  # diff颜色
  ggplot2::scale_fill_gradientn(
    colors = grDevices::colorRampPalette(c(
      "#ffffff",
      "#FFED99",
      "#85ac61",
      "#8ecde0ff",
      "#7b74e0",
      "#4941b9",
      "#991cb9"
    ))(10),
    limits = c(0, 6.5),
    breaks = seq(0, 6.5, by = 1),
    na.value = "#e9e9e9ff", # ← NA 灰色
    name = "Diff"
  ) +
  ggnewscale::new_scale("fill") +
  # 显著性颜色
  ggplot2::geom_polygon(
    ggplot2::aes(lower.x, lower.y, fill = neg_log10_p, group = group),
    colour = "white",
    linewidth = 0.1
  ) +
  ggplot2::scale_fill_gradient(
    low = "#fceeeeff",
    high = "#d65456ff",
    limits = c(0, 300),
    na.value = "#e9e9e9ff",
    name = "-log10 (P value)"
  ) +
  ggplot2::geom_text(
    data = centers,
    ggplot2::aes(x = x_center + 0.1, y = y_center - 0.12, label = label),
    size = 2.4,
    fontface = "bold"
  ) +
  ggplot2::scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    expand = c(0, 0),
    breaks = seq_along(unique(plot_df2[["comparison"]])) + 0.5,
    labels = sort(unique(plot_df2[["comparison"]])),
    sec.axis = ggplot2::dup_axis()
  ) +
  ggplot2::labs(
    title = "UCell Pos Ctrl",
    subtitle = "signif: wilcoxon rank sum test -> UCell score\n
    diff = mean score / mean score\n
    -log10 (P value) = -log10 (wilcoxon rank sum test p value)"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.y.left = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    # axis.text.x = ggplot2::element_text(vjust = 0.5, size = 10, angle = 90),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.y.left = ggplot2::element_blank(),
    axis.text.y.right = ggplot2::element_text(size = 10, face = "bold"),
    strip.text.y = ggplot2::element_text(size = 10, face = "bold"),
    strip.text.x = ggplot2::element_text(size = 8, face = "bold"),
    strip.background.y = ggplot2::element_rect(
      color = "white",
      fill = "#EEEEEE"
    ),
    strip.background.x = ggplot2::element_rect(
      color = "white",
      fill = "#EEEEEE"
    ),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
  ) +
  # 四变量分面
  ggplot2::facet_grid(
    `UCell_type` ~ pheno + sc + bulk,
    scales = "free", # must be free
    space = "free" # must be free
  )

ggplot2::ggsave(
  p,
  filename = "heatmap_combined.png",
  width = 20,
  height = 12,
  dpi = 400
)

p_hairtail <- ggplot2::ggplot(plot_df2) +
  ggplot2::geom_polygon(
    ggplot2::aes(upper.x, upper.y, fill = diff, group = group),
    colour = "grey",
    linewidth = 0.1
  ) +
  # diff颜色
  ggplot2::scale_fill_gradientn(
    colors = grDevices::colorRampPalette(c(
      "#ffffff",
      "#FFED99",
      "#85ac61",
      "#8ecde0ff",
      "#7b74e0",
      "#4941b9",
      "#991cb9"
    ))(10),
    limits = c(0, 6.5),
    breaks = seq(0, 6.5, by = 1),
    na.value = "#e9e9e9ff", # ← NA 灰色
    name = "Diff"
  ) +
  ggnewscale::new_scale("fill") +
  # 显著性颜色
  ggplot2::geom_polygon(
    ggplot2::aes(lower.x, lower.y, fill = neg_log10_p, group = group),
    colour = "white",
    linewidth = 0.1
  ) +
  ggplot2::scale_fill_gradient(
    low = "#fceeeeff",
    high = "#d65456ff",
    limits = c(0, 300),
    na.value = "#e9e9e9ff",
    name = "-log10 (P value)"
  ) +
  ggplot2::geom_text(
    data = centers,
    ggplot2::aes(x = x_center + 0.1, y = y_center - 0.12, label = label),
    size = 2.4,
    fontface = "bold"
  ) +
  ggplot2::scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    expand = c(0, 0),
    breaks = seq_along(unique(plot_df2[["comparison"]])) + 0.5,
    labels = sort(unique(plot_df2[["comparison"]])),
    sec.axis = ggplot2::dup_axis()
  ) +
  ggplot2::labs(
    title = "UCell Pos Ctrl",
    subtitle = "signif: wilcoxon rank sum test -> UCell score\n
    diff = mean score / mean score\n
    -log10 (P value) = -log10 (wilcoxon rank sum test p value)"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.y.left = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    # axis.text.x = ggplot2::element_text(vjust = 0.5, size = 10, angle = 90),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.y.left = ggplot2::element_blank(),
    axis.text.y.right = ggplot2::element_text(size = 10, face = "bold"),
    strip.text.y = ggplot2::element_text(size = 10, face = "bold"),
    strip.text.x = ggplot2::element_text(size = 8, face = "bold"),
    strip.background.y = ggplot2::element_rect(
      color = "white",
      fill = "#EEEEEE"
    ),
    strip.background.x = ggplot2::element_rect(
      color = "white",
      fill = "#EEEEEE"
    ),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
  ) +
  # 四变量分面
  ggplot2::facet_grid(
    ~ `UCell_type` + pheno + sc + bulk,
    scales = "free", # must be free
    space = "free" # must be free
  )

ggplot2::ggsave(
  p_hairtail,
  filename = "heatmap_combined_hairtail.png",
  width = 40,
  height = 8,
  dpi = 400
)
