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
  tissue = factor(tissue, levels = sort(unique(tissue))),
  x = .I
)]

data.table::setorder(spiral_f1, sample_name, method_name)
spiral_f1[, x := .I]

# 绘制
grDevices::png(
  "spiral_f1.png",
  width = 14,
  height = 10,
  units = "in",
  res = 400
)

n <- nrow(spiral_f1)

# * 5 圈
cli::cli_h2("init")

spiralize::spiral_initialize(
  xlim = c(1, n),
  start = 360,
  end = 360 * 6,
  scale_by = "curve_length",
)
spiralize::spiral_track(
  height = 0.5,
  ylim = c(0, 1),
  background_gp = grid::gpar(fill = "#EEEEEE")
)
spiralize::spiral_yaxis(side = "both", at = c(0, 0.25, 0.5, 0.75, 1))

legend_rule <- spiralize::spiral_horizon(
  spiral_f1$x,
  spiral_f1$F1,
  y_max = 1,
  use_bar = TRUE
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
sample_group_index <- spiral_f1[, .(indices = list(x)), by = sample_name]
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

palette_method_name <- c(
  "scissor" = "#B8D9A0",
  "scAB" = "#e7c0d7",
  "scPAS" = "#A8D9D8",
  "scPP" = "#F5C0A8",
  "DEGAS" = "#D5C9E8",
  "SCIPAC" = "#9BBBD9"
)
method_group_index <- spiral_f1[,
  .(indices = list(x)),
  by = .(sample_name, method_name)
]
for (i in seq_len(nrow(method_group_index))) {
  this <- method_group_index[i, ]
  cli::cli_alert_info(
    "color:{palette_method_name[this$method_name]}: {this$sample_name} - {this$method_name}"
  )
  spiralize::spiral_highlight(
    min(this$indices[[1L]]),
    max(this$indices[[1L]]),
    type = "line",
    gp = grid::gpar(col = palette_method_name[this$method_name]),
    line_width = grid::unit(6, "pt"),
  )
}

cli::cli_h2("Legend")

# * Legend
f1_lgd <- spiralize::horizon_legend(legend_rule, title = "F1 score")
method_lgd <- ComplexHeatmap::Legend(
  at = names(palette_method_name),
  title = "Method",
  legend_gp = grid::gpar(
    fill = palette_method_name
  )
)
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
  y = grid::unit(0.4, "npc"),
  just = "left"
)
ComplexHeatmap::draw(
  method_lgd,
  x = grid::unit(1, "npc") + grid::unit(2, "mm"),
  y = grid::unit(0.28, "npc"),
  just = "left"
)

grDevices::dev.off()

# rpkgkit::make_func_call_explicit(use_packages = c("ComplexHeatmap","grid","spiralize","base"))
