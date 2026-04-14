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

names(esmats) <- paste(
  basename(dirname(dirname(esmats))), # tumor type
  basename(dirname(esmats)), # bulk
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
  "SCIPAC",
  "scPP",
  "DEGAS",
  "LP_SGL",
  "PIPET"
)
esmats_data_loaded <- lapply(esmats_data, function(df) {
  cols_idx <- c(
    1,
    2,
    match(screen_methods, colnames(df))
  )
  cols_idx <- cols_idx[!is.na(cols_idx)]
  df[, cols_idx, drop = FALSE]
})

# * add sc, bulk and pheno info
esmats_data_with_info <- purrr::imap(
  esmats_data_loaded,
  function(df, name) {
    dt <- dplyr::mutate(df, agg_col = name) %>%
      tidyr::extract(
        col = agg_col,
        into = c("sc", "bulk", "pheno"),
        regex = "^([^_]+)_(.+)_(.+)$"
      ) %>%
      dplyr::rename(pos_ssGSEA = 1, neg_ssGSEA = 2) %>%
      data.table::setDT()
  },
  .progress = "ETA"
)
esmats_combined <- data.table::rbindlist(esmats_data_with_info, fill = TRUE)

esmats_combined_long <- data.table::melt(
  data = esmats_combined,
  measure.vars = screen_methods,
  variable.name = "screen_method", # 存放原列名的新列
  value.name = "screen_group", # 存放对应值的新列（可按实际含义重命名）
  na.rm = FALSE # 若需过滤缺失值，设为 TRUE
)

esmats_combined_long <- data.table::melt(
  data = esmats_combined_long,
  measure.vars = c("pos_ssGSEA", "neg_ssGSEA"),
  variable.name = "ssGSEA_type", # 存放原列名的新列
  value.name = "ssGSEA_score", # 存放对应值的新列（可按实际含义重命名）
  na.rm = FALSE # 若需过滤缺失值，设为 TRUE
)

esmats_combined_long[,
  cluster := data.table::fifelse(
    is.na(screen_group),
    NA_character_,
    cheapr::paste_(screen_method, screen_group, sep = "_")
  )
]

three_group_colors <- c(
  "Positive" = "#c24b4b",
  "Negative" = "#5189bb",
  "Neutral" = "#CECECE"
)
two_group_colors <- c(
  "Positive" = "#c24b4b",
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

# label_position <- esmats_combined_long %>%
#   group_by(cluster) %>%
#   summarise(
#     y_pos = max(ssgsea_score, na.rm = TRUE),
#     label = paste0("n=", n()),
#     .groups = "drop"
#   )

esmats_combined_long_filtered <- esmats_combined_long[
  !is.na(screen_group) & !is.na(cluster)
]
label_position <- esmats_combined_long_filtered[,
  list(
    y_pos = max(ssGSEA_score, na.rm = TRUE),
    label = paste0("n=", .N)
  ),
  by = list(cluster, sc, pheno, bulk)
]

p <- ggplot(
  esmats_combined_long_filtered,
  aes(x = cluster, y = ssGSEA_score, fill = cluster)
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
    scale = "count"
  ) +
  scale_fill_manual(values = pallete) +
  scale_y_continuous(
    breaks = scales::breaks_width(0.2),
    minor_breaks = scales::breaks_width(0.1)
  ) +
  labs(x = NULL, y = NULL) + # 分面时统一加 lab
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  ggrepel::geom_text_repel(
    data = label_position,
    aes(x = cluster, y = y_pos + 0.2, label = label),
    size = 2.8,
    fontface = "bold",
    box.padding = 0.01,
    point.padding = 0.01, # 与数据点的距离
    segment.color = NA, # 不显示连接线（如需显示可去掉这行）
    direction = "y", # 主要沿 Y 轴方向排斥
    max.overlaps = Inf # 允许尝试所有重叠情况
  ) +
  ggplot2::facet_grid(sc + bulk ~ pheno, scales = "free", space = "free")


ggplot2::ggsave(
  filename = "boxplot_combined.png",
  plot = p,
  width = 20,
  height = 24,
  dpi = 400
)
