setwd(file.path(usethis::proj_path(), "4_positive_ctrl/vis/combined"))

# library(data.table)
library(ggplot2)
# library(ggsignif)
library(patchwork)
library(gghalves)

esmats <- list.files(
  "../../esmat",
  pattern = ".qs",
  recursive = TRUE,
  full.names = TRUE
)
names(esmats) <- basename(esmats) %>%
  tools::file_path_sans_ext() %>%
  gsub("ssGSEA_score_", "", .)
names(esmats) <- paste(
  basename(dirname(dirname(esmats))),
  names(esmats),
  gsub(".*(survival|binary).*", "\\1", esmats, ignore.case = TRUE),
  sep = "_"
)

esmats_data <- lapply(esmats, function(x) {
  qs::qread(x, nthreads = 8L) # data.frame
})

# * keep only score and group
screen_methods <- c(
  "scissor",
  "scAB",
  "scPAS",
  #   "SCIPAC",
  "scPP",
  "DEGAS",
  "LP_SGL",
  "PIPET"
)
esmats_data <- lapply(esmats_data, function(df) {
  cols_idx <- c(
    1,
    2,
    match(screen_methods, colnames(df))
  )
  cols_idx <- cols_idx[!is.na(cols_idx)]
  df[, cols_idx, drop = FALSE]
})

# * add sc, bulk and pheno info
esmats_data <- purrr::imap(esmats_data, function(df, name) {
  dt <- dplyr::mutate(df, agg_col = name) %>%
    tidyr::separate(
      col = "agg_col",
      into = c("sc", "bulk", "pheno"),
      sep = "_",
      remove = FALSE
    ) %>%
    dplyr::rename(pos_ssGSEA = 1, neg_ssGSEA = 2) %>%
    data.table::setDT()
  dt[,
    (cols) := Map(function(n, v) paste0(n, " ", v), n = cols, v = .SD),
    .SDcols = screen_methods
  ]
})
esmats_combined <- data.table::rbindlist(esmats_data)
esmats_combined_long <- data.table::melt(
  data = esmats_combined,
  measure.vars = screen_methods,
  variable.name = "screen_method", # 存放原列名的新列
  value.name = "screen_group", # 存放对应值的新列（可按实际含义重命名）
  na.rm = FALSE # 若需过滤缺失值，设为 TRUE
)

three_group_colors <- c(
  "Positive" = "#ff3333",
  "Negative" = "#386c9b",
  "Neutral" = "#CECECE"
)
two_group_colors <- c(
  "Positive" = "#ff3333",
  "Other" = "#CECECE"
)

pallete <- c(
  rep(two_group_colors, 2),
  rep(three_group_colors, 6)
)
names(pallete) <- c(
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

p <- ggplot(
  esmats_combined_long,
  aes(x = cluster, y = ssgsea_score, fill = cluster)
) +
  geom_half_boxplot(
    outlier.alpha = 0.2,
    outlier.size = 0.5,
    outlier.colour = "#cececeff",
    width = 0.65,
    alpha = 0.8,
    side = "l",
    show.legend = FALSE
  ) +
  geom_half_violin(
    side = "r",
    trim = FALSE,
    alpha = 0.5,
    show.legend = FALSE,
    width = 0.65,
    scale = "width"
  ) +
  scale_fill_manual(values = pallete, alpha = 0.5) +
  scale_y_continuous(
    breaks = scales::breaks_width(0.2),
    minor_breaks = scales::breaks_width(0.1)
  ) +
  labs(x = NULL, y = NULL) + # 分面时统一加 lab
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  ggrepel::geom_text_repel(
    data = label_position,
    aes(x = cluster, y = y_pos + 0.05, label = label),
    size = 3,
    fontface = "bold",
    box.padding = 0.5, # 文字周围的填充空间
    point.padding = 0.5, # 与数据点的距离
    segment.color = NA, # 不显示连接线（如需显示可去掉这行）
    direction = "y", # 主要沿 Y 轴方向排斥
    max.overlaps = Inf # 允许尝试所有重叠情况
  ) +
  ggplot2::facet_grid(sc + bulk ~ pheno, scales = "free", space = "free")


ggplot2::ggsave(
  filename = "boxplot_combined.png",
  plot = p,
  width = 10,
  height = 10,
  dpi = 400
)
