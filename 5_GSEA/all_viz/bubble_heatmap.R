# ! run after `combine_data.R`
library(data.table)

setwd(file.path(usethis::proj_path(), "5_GSEA/all_viz"))

gsea_score_combined <- data.table::fread(
  "/home/data/sigbridger/GSEA/mean_score_gsea/gsea_score_combined.csv"
)

gsea_score_list <- gsea_score_combined[, .(list(.SD)), by = data_name]$V1
names(gsea_score_list) <- unique(gsea_score_combined$data_name)

# ============================================================
# Bubble heatmap: phenotype, sc, bulk 分面
# ============================================================

future::plan(future.mirai::mirai_multisession(), workers = 4L)
dir.create("plots", recursive = TRUE, showWarnings = FALSE)


furrr::future_walk(
  .x = gsea_score_list,
  .f = \(dt) {
    # * 熔化为长格式：每个基因集作为单独行
    id_vars <- c(
      "method",
      "status",
      "count",
      "phenotype",
      "sc",
      "bulk"
    )
    value_cols <- setdiff(colnames(dt), id_vars)

    gsea_long <- data.table::melt(
      dt,
      id.vars = id_vars,
      measure.vars = value_cols,
      variable.name = "geneset",
      value.name = "score"
    )

    # * 缩写基因集名称（便于显示）
    # 1. 保留方向性修饰词，用箭头/缩写替代
    gsea_long[, geneset_short := geneset]
    gsea_long[,
      geneset_short := gsub(
        "Genes up-regulated",
        "UP:",
        geneset_short,
        ignore.case = TRUE
      )
    ]
    gsea_long[,
      geneset_short := gsub(
        "Genes down-regulated",
        "DOWN:",
        geneset_short,
        ignore.case = TRUE
      )
    ]
    gsea_long[,
      geneset_short := gsub(
        "up-regulated",
        "UP:",
        geneset_short,
        ignore.case = TRUE
      )
    ]
    gsea_long[,
      geneset_short := gsub(
        "down-regulated",
        "DOWN:",
        geneset_short,
        ignore.case = TRUE
      )
    ]
    # 2. 去除 "Genes " 前缀（普通描述性前缀）
    gsea_long[,
      geneset_short := gsub(
        "Genes (encoding |involved in |defining |important for |regulating |mediating )",
        "",
        geneset_short,
        ignore.case = TRUE
      )
    ]
    gsea_long[,
      geneset_short := gsub(
        "Genes specifically ",
        "",
        geneset_short,
        ignore.case = TRUE
      )
    ]
    gsea_long[,
      geneset_short := gsub(
        "^Genes ",
        "",
        geneset_short
      )
    ]
    # 3. 去除其他冗余前缀
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
        "^genes ",
        "",
        geneset_short
      )
    ]
    # 4. 去除长注释后缀
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

    # * method+status 组合标签 + count 映射
    gsea_long[, method_status := paste0(method, "_", status)]
    count_lab <- gsea_long[, .(count = unique(count)), by = method_status]
    y_labels <- setNames(
      paste0(count_lab$method_status, "\n(n = ", count_lab$count, ")"),
      count_lab$method_status
    )

    sc <- unique(dt[, sc]) # 单细胞
    bulk <- unique(dt[, bulk]) #  bulk
    phenotype <- unique(dt[, phenotype]) #  表型

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
      ggplot2::geom_point() +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", score)),
        size = 3,
        color = "grey10",
        vjust = -1.5,
        show.legend = FALSE
      ) +
      ggplot2::scale_y_discrete(labels = y_labels) +
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
        range = c(1, 6),
        name = "|Score|",
        limits = c(0, 1)
      ) +
      ggplot2::labs(
        x = "Hallmark Gene Set",
        y = "Method × Status",
        title = glue::glue("{phenotype} {sc} {bulk}")
      ) +
      cowplot::theme_cowplot() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = ggplot2::element_text(size = 14),
        panel.grid.major = ggplot2::element_line(
          color = "grey90",
          linewidth = 0.3
        ),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.key.size = ggplot2::unit(0.8, "lines")
      ) +
      ggplot2::guides(
        color = ggplot2::guide_colorbar(
          order = 1,
          barwidth = 12,
          barheight = 0.6
        ),
        size = ggplot2::guide_legend(order = 2, nrow = 1)
      )

    # * 保存
    plot_name <- glue::glue(
      "{phenotype}_{sc}_{bulk}_bubble_hallmark_score.png"
    )
    ggplot2::ggsave(
      file.path("plots", plot_name),
      p,
      width = 18,
      height = 14,
      limitsize = FALSE
    )
  },
  .progress = TRUE
)
