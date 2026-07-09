# ! run this after `mean_score.R`
# ! all tissue
library(dplyr)
library(data.table)

setwd(file.path(usethis::proj_path(), "5_GSEA/all_viz"))

data_dir <- "/home/data/sigbridger/GSEA/mean_score_gsea"
files <- list.files(
  data_dir,
  full.names = TRUE,
  include.dirs = FALSE
)
names(files) <- basename(files) %>%
  tools::file_path_sans_ext() %>%
  gsub("_mean.*", "", .)

gsea_score_loaded <- purrr::imap(
  files,
  ~ {
    cli::cli_alert_info("{.file {.y}} loaded")
    loaded <- data.table::fread(.x)
    loaded[, data_name := .y]
  }
)
gsea_score_combinedd <- dplyr::bind_rows(gsea_score_loaded)
data.table::set(
  gsea_score_combinedd,
  j = c("phenotype"),
  value = gsub("_.*", "", gsea_score_combinedd[, data_name])
)
data.table::set(
  gsea_score_combinedd,
  j = c("sc"),
  value = gsub(
    ".*(her2|tnbc|ov|lung).*",
    "\\1",
    gsea_score_combinedd[, data_name],
    ignore.case = TRUE
  )
)
data.table::set(
  gsea_score_combinedd,
  j = c("bulk"),
  value = gsub("^([^_]*_){2}", "", gsea_score_combinedd[, data_name])
)

# ============================================================
# Bubble heatmap: phenotype, sc, bulk 分面
# ============================================================

# * 熔化为长格式：每个基因集作为单独行
id_vars <- c(
  "method",
  "status",
  "count",
  "phenotype",
  "sc",
  "bulk",
  "data_name"
)
value_cols <- setdiff(colnames(gsea_score_combinedd), id_vars)

gsea_long <- data.table::melt(
  gsea_score_combinedd,
  id.vars = id_vars,
  measure.vars = value_cols,
  variable.name = "geneset",
  value.name = "score"
)

# * 缩写基因集名称（便于显示）
gsea_long[,
  geneset_short := gsub(
    "Genes (encoding |involved in |defining |up-regulated |down-regulated |important for |regulating |specific|up regulated|down regulated)",
    "",
    geneset,
    ignore.case = TRUE
  )
]
gsea_long[,
  geneset_short := gsub(
    "(A subgroup of |processing of |production of |formation of |development of |metabolism of |components of )",
    "",
    geneset_short,
    ignore.case = TRUE
  )
]
gsea_long[,
  geneset_short := gsub(
    "(, as in.*|, e\\.g\\..*|, a cellular.*|, which.*|; also.*)",
    "",
    geneset_short
  )
]
# 截断过长的名称
gsea_long[,
  geneset_short := ifelse(
    nchar(geneset_short) > 45,
    paste0(substr(geneset_short, 1, 42), "..."),
    geneset_short
  )
]

# * method+status 组合标签
gsea_long[, method_status := paste0(method, "_", status)]

# * 气泡热图
p <- ggplot2::ggplot(
  gsea_long,
  ggplot2::aes(
    x = geneset_short,
    y = method_status,
    color = score,
    size = abs(score)
  )
) +
  ggplot2::geom_point(alpha = 0.85) +
  ggplot2::geom_text(
    ggplot2::aes(label = count),
    size = 2.2,
    color = "grey30",
    vjust = -1.6,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Mean Score"
  ) +
  ggplot2::scale_size_continuous(
    range = c(1, 3),
    name = "|Score|",
    limits = c(0, 1)
  ) +
  ggplot2::facet_grid(
    phenotype + bulk ~ sc,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::labs(x = "Hallmark Gene Set", y = "Method × Status") +
  cowplot::theme_cowplot() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = ggplot2::element_text(size = 8),
    strip.text = ggplot2::element_text(size = 7),
    strip.background = ggplot2::element_rect(fill = "grey92", color = NA),
    panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.size = ggplot2::unit(0.8, "lines")
  ) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(order = 1, barwidth = 12, barheight = 0.6),
    size = ggplot2::guide_legend(order = 2, nrow = 1)
  )

# * 保存
dir.create("plots", recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
  file.path("plots", "bubble_hallmark_gsea_score.png"),
  p,
  width = 22,
  height = 18,
  limitsize = FALSE
)
cli::cli_alert_success("Bubble heatmap saved.")
