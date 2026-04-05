setwd(file.path(usethis::proj_path(), "4_positive_ctrl/vis/combined"))

library(data.table)

data_survival <- list.files(
  path = file.path(usethis::proj_path(), "4_positive_ctrl/vis/survival"),
  pattern = "*.csv",
  full.names = TRUE,
  recursive = TRUE
)
data_binary <- list.files(
  path = file.path(usethis::proj_path(), "4_positive_ctrl/vis/binary"),
  pattern = "*.csv",
  full.names = TRUE,
  recursive = TRUE
)
data <- c(data_survival, data_binary)
names(data) <- basename(data) %>% tools::file_path_sans_ext()

data_loaded <- lapply(X = data, FUN = data.table::fread)
data_loaded <- purrr::imap(.x = data_loaded, .f = \(dt, name) {
  dt[,
    tumor_type := stringi::stri_match_first_regex(
      name,
      "^[^_]*_([^_]*)_"
    )[, 2] %>%
      toupper()
  ]
})
data_combined <- data.table::rbindlist(data_loaded)
data_combined <- data_combined[group.y == "Positive vs Rest"]

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
plot_df <- data_combined %>%
  tidyr::unite(
    "comparison",
    screen_method,
    group.y,
    remove = FALSE
  ) %>%
  dplyr::right_join(method_comparisons, by = c("bulk", "comparison")) %>%
  dplyr::mutate(label = ifelse(is.na(label), "NA", label))

# 在调用 triangle_data 之前
plot_df <- plot_df %>%
  dplyr::mutate(
    # 将 bulk 和 comparison 转换为因子，并指定你想要的顺序
    comparison = factor(
      comparison,
      levels = unique(sort(plot_df$comparison))
    ),
    bulk = factor(bulk, levels = unique(sort(plot_df$bulk)))
  ) %>%
  dplyr::arrange(bulk, comparison)

plot_df <- dplyr::filter(
  plot_df,
  !is.na(tumor_type) & !is.na(type_pheno) & !is.na(bulk)
)

plot_df2 <- triangle_data(plot_df, row = "bulk", col = "comparison")

plot_df2 <- plot_df2 %>%
  dplyr::filter(
    !is.na(tumor_type) & !is.na(type_pheno) & !is.na(bulk)
  )

# # ? significance label
# points <- plot_df2
# points$x <- rep(0, nrow(points))
# points$y <- rep(0, nrow(points))

# ? significance label position
centers <- plot_df2 %>%
  dplyr::group_by(
    bulk,
    comparison,
    type_pheno,
    `ssGSEA type`,
    tumor_type
  ) %>%
  dplyr::summarise(
    x_center = mean(lower.x),
    y_center = mean(lower.y),
    label = dplyr::first(label),
    .groups = "drop"
  ) %>%
  dplyr::filter(label != "" & !is.na(label))

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
      "#3c3597",
      "#465d9bff"
    ))(10),
    limits = c(0.5, 3),
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
    breaks = seq_along(unique(plot_df2[["bulk"]])) + 0.5,
    expand = c(0, 0),
    labels = sort(unique(plot_df2[["bulk"]]))
  ) +
  ggplot2::scale_y_continuous(
    expand = c(0, 0),
    breaks = seq_along(unique(plot_df2[["comparison"]])) + 0.5,
    labels = sort(unique(plot_df2[["comparison"]])),
    sec.axis = ggplot2::dup_axis()
  ) +
  ggplot2::labs(
    title = "ssGSEA Pos Ctrl",
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.y.left = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    # axis.text.x = ggplot2::element_text(vjust = 0.5, size = 10, angle = 90),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.y.left = ggplot2::element_blank(),
    axis.text.y.right = ggplot2::element_text(size = 10, face = "bold"),
    strip.text.y = ggplot2::element_text(size = 7),
    strip.text.x = ggplot2::element_text(size = 7),
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
    type_pheno + `ssGSEA type` ~ tumor_type + bulk,
    scales = "free",
    space = "free"
  )

ggplot2::ggsave(
  p,
  filename = "heatmap_combined.png",
  width = 10,
  height = 10,
  dpi = 4000
)
