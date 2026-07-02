setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ============================
# 螺旋图：6 个 method 各占 60°，tissue 在螺旋轴上，F1 为径向值
# ============================

combined2 <- data.table::fread("combined.csv")


# 按 method_name + tissue 汇总 F1 均值
spiral_f1 <- combined2[, c(
  "screen_method",
  "tissue",
  "F1",
  "method_name",
  "arg_sample"
)]

# 排序：method → tissue，生成连续 x 索引

spiral_f1[, `:=`(
  sample_id = gsub(".*samples(.).*", "\\1", screen_method)
)][,
  `:=`(
    sample_id = data.table::fifelse(
      stringr::str_detect(sample_id, "sample"),
      "1",
      sample_id
    )
  )
][, `:=`(
  sample_name = paste(tissue, "dataset", sample_id),
  method_name = factor(method_name, levels = sort(unique(method_name))),
  tissue = factor(tissue, levels = sort(unique(tissue)))
)]

spiral_f1_by_method <- split(spiral_f1, by = "method_name")

spiral_f1_by_method <- lapply(spiral_f1_by_method, \(dt) {
  data.table::setorder(dt, sample_name, method_name)
  dt[, x := .I]
})

purrr::iwalk(
  spiral_f1_by_method,
  \(dt, dt_name) {
    # 绘制
    grDevices::png(
      glue::glue("spiral_f1_{dt_name}.png"),
      width = 14,
      height = 10,
      units = "in",
      res = 400
    )

    n <- nrow(dt)

    # * 5 圈
    cli::cli_h2("init")

    spiralize::spiral_initialize(
      xlim = c(1, n),
      start = 360,
      end = 360 * 5,
      scale_by = "curve_length",
    )
    spiralize::spiral_track(
      height = 0.5,
      ylim = c(0, 1),
      background_gp = grid::gpar(fill = "#EEEEEE")
    )
    spiralize::spiral_yaxis(side = "both", at = c(0, 0.25, 0.5, 0.75, 1))

    col_fun <- circlize::colorRamp2(c(0, 1), c("#ffdddd", "red"))
    spiralize::spiral_rect(
      seq_len(n - 1),
      0,
      2:n,
      1,
      gp = grid::gpar(fill = col_fun(dt$F1), col = NA)
    )

    grid::grid.rect(gp = grid::gpar(fill = NA))
    grid::grid.text(
      "F1 score",
      0,
      0,
      default.units = "native",
      gp = grid::gpar(fontfamily = "bold")
    )

    cli::cli_h2("Annotation")

    # * Group Annotation
    sample_group_index <- dt[, .(indices = list(x)), by = sample_name]
    palette_sample_name <- c(
      "brca_her2 dataset 1" = "#00b8c2",
      "brca_her2 dataset 2" = "#0067a3",
      "brca_tnbc dataset 1" = "#009940",
      "brca_tnbc dataset 2" = "#9b008e",
      "lung dataset 1" = "#9b0008"
    )
    for (i in seq_len(nrow(sample_group_index))) {
      this <- sample_group_index[i, ]
      cli::cli_alert_info(
        "color:{palette_sample_name[this$sample_name]}: {this$sample_name}"
      )
      spiralize::spiral_highlight(
        min(this$indices[[1L]]),
        max(this$indices[[1L]]),
        type = "line",
        gp = grid::gpar(col = palette_sample_name[this$sample_name]),
        line_width = grid::unit(12, "pt")
      )
    }

    cli::cli_h2("Legend")

    # * Legend
    f1_lgd <- ComplexHeatmap::Legend(title = "F1 score", col_fun = col_fun)
    sample_lgd <- ComplexHeatmap::Legend(
      at = names(palette_sample_name),
      title = "Sample",
      legend_gp = grid::gpar(
        fill = palette_sample_name
      )
    )

    cli::cli_h2("Draw")

    ComplexHeatmap::draw(
      f1_lgd,
      x = grid::unit(1, "npc") + grid::unit(2, "mm"),
      just = "left"
    )
    ComplexHeatmap::draw(
      sample_lgd,
      x = grid::unit(1, "npc") + grid::unit(2, "mm"),
      y = grid::unit(0.35, "npc"),
      just = "left"
    )

    grDevices::dev.off()
  },
  .progress = "Drawing"
)

# rpkgkit::make_func_call_explicit(use_packages = c("ComplexHeatmap","grid","spiralize","base"))
