setwd(file.path(usethis::proj_path(), "2_method_acc/viz"))

library(dplyr)
library(data.table)

arg_sample_dirs <- c(
  "../brca_her2",
  "../brca_tnbc",
  "../lung"
)

arg_samples <- purrr::map(
  arg_sample_dirs,
  ~ {
    l <- list.files(
      .x,
      pattern = ".*arg_samples.*",
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )
    names(l) <- basename(l) %>% tools::file_path_sans_ext()
    return(l)
  }
)
names(arg_samples) <- basename(arg_sample_dirs)


loaded_arg_samples <- purrr::map(
  arg_samples,
  ~ purrr::map(.x, data.table::fread)
)

metrics_cols <- c(
  "TPR",
  "FPR",
  "TNR",
  "FNR",
  "Precision",
  "Recall",
  "F1",
  "Accuracy"
)

concreted_arg_samples <- purrr::map(
  loaded_arg_samples,
  ~ purrr::map(.x, function(dt) {
    arg_cols <- dplyr::setdiff(colnames(dt), metrics_cols)

    for (col in arg_cols) {
      data.table::set(dt, j = col, value = paste0(col, "=", dt[[col]]))
    }
    dt[, arg_sample := do.call(paste, c(.SD, sep = "| ")), .SDcols = arg_cols]
    dt[, (arg_cols) := NULL]
    dt
  })
)

combined <- purrr::imap(
  concreted_arg_samples,
  ~ {
    dt_tissue <- purrr::imap(.x, function(dt, name) {
      dt[, screen_method := name]
      return(dt)
    }) %>%
      dplyr::bind_rows()

    data.table::set(dt_tissue, j = "tissue", value = .y)
    dt_tissue
  }
) %>%
  dplyr::bind_rows()

combined2 <- combined %>%
  dplyr::mutate(
    # 判断条件：是否包含LP_SGL
    method_name = stringr::str_extract(
      screen_method,
      "scissor|scab|scpas|scpp|degas|lp_sgl|scipac|pipet"
    ),
    method_name = dplyr::case_when(
      method_name == "scab" ~ "scAB",
      method_name == "scpas" ~ "scPAS",
      method_name == "scpp" ~ "scPP",
      method_name == "degas" ~ "DEGAS",
      method_name == "lp_sgl" ~ "LP_SGL",
      method_name == "scipac" ~ "SCIPAC",
      method_name == "pipet" ~ "PIPET",
      TRUE ~ method_name
    )
  )

data.table::fwrite(combined2, "combined.csv")


# ----------------------------------------------------------------------------------------------------------
# * Viz
combined2 <- data.table::fread("combined.csv")

# 定义完整方法顺序
all_methods <- c(
  "scissor",
  "scAB",
  "scPAS",
  "scPP",
  "DEGAS",
  "LP_SGL",
  "SCIPAC",
  "PIPET"
)

# 为缺失的方法（LP_SGL、PIPET 等可能无数据）添加 NA 占位行
existing_methods <- unique(combined2$method_name)
missing_methods <- setdiff(all_methods, existing_methods)
if (length(missing_methods) > 0) {
  placeholder <- combined2[rep(1, length(missing_methods))]
  for (j in seq_along(placeholder)) {
    data.table::set(
      placeholder,
      j = j,
      value = rep(NA_real_, length(missing_methods))
    )
  }
  placeholder[, method_name := missing_methods]
  combined2 <- data.table::rbindlist(list(combined2, placeholder), fill = TRUE)
}

combined2[, method_name := factor(method_name, levels = all_methods)]

palette <- c(
  "scissor" = "#B8D9A0",
  "scAB" = "#E8E0DB",
  "scPAS" = "#A8D9D8",
  "scPP" = "#F5C0A8",
  "DEGAS" = "#D5C9E8",
  "LP_SGL" = "#F9D586",
  "SCIPAC" = "#9BBBD9",
  "PIPET" = "#E2A9B8"
)

p_f1 <- combined2 %>%
  ggplot2::ggplot(ggplot2::aes(x = method_name, y = F1)) +
  # 极值线：每组最小到最大的范围，带横杠帽（空组返回 NA 避免报错）
  ggplot2::stat_summary(
    fun.data = function(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0) {
        return(data.frame(ymin = NA_real_, ymax = NA_real_, y = NA_real_))
      }
      data.frame(ymin = min(x), ymax = max(x), y = mean(x))
    },
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.8,
    color = "grey30"
  ) +
  # 箱线图，按 method_name 填充 palette 颜色，中位数线调细
  ggplot2::geom_boxplot(
    ggplot2::aes(fill = method_name),
    alpha = 0.5,
    width = 0.6,
    outlier.shape = NA,
    fatten = 1,
    na.rm = TRUE
  ) +
  # 散点按 method_name 分组着色
  ggplot2::geom_jitter(
    ggplot2::aes(color = method_name),
    size = 1.5,
    alpha = 0.4,
    width = 0.2
  ) +
  # 应用调色板
  ggplot2::scale_fill_manual(values = palette, guide = "none") +
  ggplot2::scale_color_manual(
    values = c(
      "scissor" = "#6a914d",
      "scAB" = "#9f765d",
      "scPAS" = "#3d8887",
      "scPP" = "#9c644b",
      "DEGAS" = "#513b74",
      "LP_SGL" = "#7d642d",
      "SCIPAC" = "#2a547d",
      "PIPET" = "#732b3e"
    ),
    guide = "none"
  ) +
  ggplot2::scale_x_discrete(drop = FALSE) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  cowplot::theme_cowplot() +
  ggplot2::labs(
    title = "Scatter-Boxplot of F1 by Method and Tissue",
    x = "",
    y = "F1"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 12),
    axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
    axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      color = "#EEEEEE",
      linewidth = 0.1
    ),
    strip.text = ggplot2::element_text(face = "bold", size = 11),
    strip.background = ggplot2::element_rect(
      fill = "#EEEEEE",
      color = "white"
    )
  )

ggplot2::ggsave("p_f1.png", p_f1, width = 8, height = 6, dpi = 400)

p_acc <- combined2 %>%
  ggplot2::ggplot(ggplot2::aes(x = method_name, y = Accuracy)) +
  # 极值线：每组最小到最大的范围，带横杠帽（空组返回 NA 避免报错）
  ggplot2::stat_summary(
    fun.data = function(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0) {
        return(data.frame(ymin = NA_real_, ymax = NA_real_, y = NA_real_))
      }
      data.frame(ymin = min(x), ymax = max(x), y = mean(x))
    },
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.8,
    color = "grey30"
  ) +
  # 箱线图，按 method_name 填充 palette 颜色，中位数线调细
  ggplot2::geom_boxplot(
    ggplot2::aes(fill = method_name),
    alpha = 0.5,
    width = 0.6,
    outlier.shape = NA,
    median.linewidth = 1,
    na.rm = TRUE
  ) +
  # 散点按 method_name 分组着色
  ggplot2::geom_jitter(
    ggplot2::aes(color = method_name),
    size = 1.5,
    alpha = 0.4,
    width = 0.2
  ) +
  # 应用调色板
  ggplot2::scale_fill_manual(values = palette, guide = "none") +
  ggplot2::scale_color_manual(
    values = c(
      "scissor" = "#6a914d",
      "scAB" = "#9f765d",
      "scPAS" = "#3d8887",
      "scPP" = "#9c644b",
      "DEGAS" = "#513b74",
      "LP_SGL" = "#7d642d",
      "SCIPAC" = "#2a547d",
      "PIPET" = "#732b3e"
    ),
    guide = "none"
  ) +
  ggplot2::scale_x_discrete(drop = FALSE) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  cowplot::theme_cowplot() +
  ggplot2::labs(
    title = "Scatter-Boxplot of Accuracy by Method and Tissue",
    x = "",
    y = "Accuracy"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 12),
    axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
    axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      color = "#EEEEEE",
      linewidth = 0.1
    ),
    strip.text = ggplot2::element_text(face = "bold", size = 11),
    strip.background = ggplot2::element_rect(
      fill = "#EEEEEE",
      color = "white"
    )
  )

ggplot2::ggsave("p_acc.png", p_acc, width = 8, height = 6, dpi = 400)
