setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(bench) # ! This is necessary
library(ggplot2)
library(dplyr)

output_dir <- "/home/data/sigbridger/method_time_cost"
screen_methods <- list.files(output_dir) %>% grepv("\\.qs", ., invert = TRUE)

# ? read files
purrr::walk(
  screen_methods,
  function(method) {
    b_path <- file.path(output_dir, method, paste0(method, "_binary_bench.qs"))
    s_path <- file.path(
      output_dir,
      method,
      paste0(method, "_survival_bench.qs")
    )

    if (!file.exists(b_path)) {
      cli::cli_warn(
        "(binary) Files not found for method: {method}, skipping..."
      )
      return(NULL)
    }

    b <- qs::qread(b_path, nthreads = 2L)
    # add a column to identify method
    b$method <- method

    assign(
      paste0(method, "_binary_tweaked"),
      b,
      envir = .GlobalEnv
    )

    if (!file.exists(s_path)) {
      cli::cli_warn(
        "(survivl) Files not found for method: {method}, skipping..."
      )
      return(NULL)
    }

    s <- qs::qread(s_path, nthreads = 2L)
    s$method <- method
    assign(
      paste0(method, "_survival_tweaked"),
      s,
      envir = .GlobalEnv
    )

    gc(verbose = FALSE)
  },
  .progress = "reading"
)

# ------ combined data ---------------

survival_tweaked <- dplyr::bind_rows(
  purrr::map(
    screen_methods,
    function(method) {
      df <- get0(paste0(method, "_survival_tweaked"), envir = .GlobalEnv)
      df$expression <- NULL
      df
    }
  )
)

binary_tweaked <- dplyr::bind_rows(
  purrr::map(
    screen_methods,
    function(method) {
      df <- get0(paste0(method, "_binary_tweaked"), envir = .GlobalEnv)
      df$expression <- NULL
      df
    }
  )
)

rm(list = ls(pattern = "_survival|_binary"), envir = .GlobalEnv)

qs::qsave(
  survival_tweaked,
  file.path(output_dir, "combined_survival_bench.qs"),
  nthreads = 64L
)

qs::qsave(
  binary_tweaked,
  file.path(output_dir, "combined_binary_bench.qs"),
  nthreads = 64L
)

# survival_tweaked <- qs::qread(
#   file.path(output_dir, "combined_survival_bench.qs"),
#   nthreads = 2L
# )
# binary_tweaked <- qs::qread(
#   file.path(output_dir, "combined_binary_bench.qs"),
#   nthreads = 2L
# )

# ---- survival plot ----

surv_plot_df <- survival_tweaked %>%
  tidyr::unnest(c(time, gc))
# %>%  dplyr::mutate(expression = as.character(expression))

# 计算每组的均值并保留用于标注的 n_cells（只取第一条），并按 method, mem_alloc 排序
group_stats <- surv_plot_df %>%
  dplyr::group_by(method, n_cells) %>%
  dplyr::summarise(
    n_cells = as.character(dplyr::first(n_cells)),
    mem_alloc = mean(as.numeric(mem_alloc), na.rm = TRUE),
    time = mean(as.numeric(time), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(method, mem_alloc)

surv_p <- surv_plot_df %>%
  ggplot2::ggplot(ggplot2::aes(
    x = mem_alloc,
    y = time,
    color = method # 👈 改为 method 为主色变量
  )) +
  # 原始数据点（可保留或淡化；如只关心均值，可注释掉）
  ggplot2::geom_point(alpha = 0.3, size = 1) +

  # 每个 method 的均值连线 & 均值点（带白心）
  ggplot2::geom_line(
    data = group_stats,
    ggplot2::aes(group = method), # group 仍按 method 分组
    linewidth = 1
  ) +
  ggplot2::geom_point(
    data = group_stats,
    ggplot2::aes(),
    shape = 21,
    fill = "white",
    color = "black",
    size = 3
  ) +

  # 标注（n_cells），避免重叠
  ggrepel::geom_text_repel(
    data = group_stats,
    ggplot2::aes(
      label = paste0(as.numeric(n_cells) / 1000, "k")
    ),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +

  bench::scale_color_bench_expr(scales::brewer_pal(
    type = "qual",
    palette = 3
  )) +

  ggplot2::labs(
    title = "Screening Analysis Benchmark",
    subtitle = "Sampling sc GSE165897; bulk GSE140082; Survival as phenotype",
    x = "Memory Allocation",
    y = "Time",
    color = "Method"
  ) +

  # 移除 legend.position = "none" 或设为 default
  cowplot::theme_cowplot() +
  ggplot2::theme(
    legend.position = "right", # 或 "top"/"bottom"
    legend.title = ggplot2::element_text(face = "bold")
  )

# surv_p

# 保存并返回绘图对象
ggplot2::ggsave(
  file.path("surv_bench.png"),
  plot = surv_p,
  width = 7,
  height = 5,
  dpi = 300
)

# ---- binary plot ----

binary_plot_df <- binary_tweaked %>%
  tidyr::unnest(c(time, gc))
# %>%  dplyr::mutate(expression = as.character(expression))

# 计算每组的均值并保留用于标注的 n_cells（只取第一条），并按 method, mem_alloc 排序
binary_group_stats <- binary_plot_df %>%
  dplyr::group_by(method, mem_alloc) %>%
  dplyr::summarise(
    n_cells = as.character(dplyr::first(n_cells)),
    mem_alloc = mean(as.numeric(mem_alloc), na.rm = TRUE),
    time = mean(as.numeric(time), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(method, mem_alloc)

bin_p <- binary_plot_df %>%
  ggplot2::ggplot(ggplot2::aes(
    x = mem_alloc,
    y = time,
    color = method # 👈 改为 method 为主色变量
  )) +
  # 原始数据点（可保留或淡化；如只关心均值，可注释掉）
  ggplot2::geom_point(alpha = 0.3, size = 1) +

  # 每个 method 的均值连线 & 均值点（带白心）
  ggplot2::geom_line(
    data = binary_group_stats,
    ggplot2::aes(group = method), # group 仍按 method 分组
    linewidth = 1
  ) +
  ggplot2::geom_point(
    data = binary_group_stats,
    ggplot2::aes(),
    shape = 21,
    fill = "white",
    color = "black",
    size = 3
  ) +

  # 标注（n_cells），避免重叠
  ggrepel::geom_text_repel(
    data = binary_group_stats,
    ggplot2::aes(
      label = paste0(as.numeric(n_cells) / 1000, "k")
    ),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +

  bench::scale_color_bench_expr(scales::brewer_pal(
    type = "qual",
    palette = 3
  )) +

  ggplot2::labs(
    title = "Screening Analysis Benchmark",
    subtitle = "Sampling sc GSE165897; bulk GSE140082; Binary phenotype",
    x = "Memory Allocation",
    y = "Time",
    color = "Method"
  ) +

  # 移除 legend.position = "none" 或设为 default
  cowplot::theme_cowplot() +
  ggplot2::theme(
    legend.position = "right", # 或 "top"/"bottom"
    legend.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  file.path("binary_bench.png"),
  plot = bin_p,
  width = 7,
  height = 5,
  dpi = 300
)

cli::cli_h1("All done!")
