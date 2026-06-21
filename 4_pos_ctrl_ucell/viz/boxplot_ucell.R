setwd(file.path(usethis::proj_path(), "4_pos_ctrl_ucell/viz"))

library(data.table)
library(ggplot2)
# library(ggsignif)
library(patchwork)
library(gghalves)

esmats <- qs::qread("../diff_test/method_labels.qs", nthreads = 2L)
matched_data <- qs::qread("../diff_test/matched_data.qs", nthreads = 2L)

# ? Find matched datasets
find_name <- function(chr = character) {
  sc <- gsub(
    ".*(her2|tnbc|lung|ov).*",
    "\\1",
    chr,
    ignore.case = TRUE
  )
  bulk <- tolower(gsub(
    ".*(TCGA.*|GSE.*)",
    "\\1",
    chr,
    ignore.case = TRUE
  ))
  pheno <- gsub(
    "(survival|binary).*",
    "\\1",
    chr,
    ignore.case = TRUE
  )
  list(
    sc = sc,
    bulk = bulk,
    pheno = pheno
  )
}


combined_data <- purrr::map2(
  matched_data,
  esmats,
  function(x, y) {
    dplyr::bind_cols(x[1:2], y)
  },
  .progress = "Combining"
)

combined_long <- purrr::imap(
  combined_data,
  function(x, y) {
    info <- find_name(y)
    dt <- data.table::setDT(x)
    col_names <- colnames(dt)
    n_cols <- ncol(dt)

    dt_long <- data.table::melt(
      data = dt,
      measure.vars = col_names[3:n_cols],
      variable.name = "screen_method", # 存放原列名的新列
      value.name = "screen_group", # 存放对应值的新列（可按实际含义重命名）
      na.rm = FALSE # 若需过滤缺失值，设为 TRUE
    )

    dt_long <- data.table::melt(
      data = dt_long,
      measure.vars = col_names[1:2],
      variable.name = "UCell_type", # 存放原列名的新列
      value.name = "UCell_score", # 存放对应值的新列（可按实际含义重命名）
      na.rm = FALSE # 若需过滤缺失值，设为 TRUE
    )
    dt_long[,
      `:=`(
        cluster = data.table::fifelse(
          is.na(screen_group),
          NA_character_,
          cheapr::paste_(screen_method, screen_group, sep = "_")
        ),
        sc = info$sc,
        bulk = info$bulk,
        pheno = info$pheno
      )
    ][, count := .N / 2, by = cluster]
    dt_long
  },
  .progress = "Melting"
)

# * contains NA
all_combined <- data.table::rbindlist(combined_long) # ~ 10M rows

# NA_dt <- all_combined[is.na(screen_group), ]

label_position <- all_combined[,
  list(
    y_pos = max(UCell_score, na.rm = TRUE),
    label = paste0("n=", .N / 2)
  ),
  by = list(cluster, sc, pheno, bulk)
]

qs::qsave(all_combined, file = "all_combined.qs", nthreads = 8L) # < 40MB
# all_combined <- qs::qread("all_combined.qs", nthreads = 8L)

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


p <- ggplot(
  all_combined,
  aes(x = cluster, y = UCell_score, fill = cluster)
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
    strip.text = element_text(face = "bold", size = 11),
    strip.background.y = ggplot2::element_rect(
      color = "white",
      fill = "#EEEEEE"
    ),
    strip.background.x = ggplot2::element_rect(
      color = "white",
      fill = "#EEEEEE"
    ),
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
  ggplot2::facet_grid(sc + bulk ~ pheno, scales = "free", space = "fixed")


ggplot2::ggsave(
  filename = "boxplot_combined.png",
  plot = p,
  width = 20,
  height = 24,
  dpi = 400
)

cli::cli_h1("Boxplot done")
gc()
