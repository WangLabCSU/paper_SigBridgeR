setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(data.table)

esmat_root <- "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/esmat/binary/ov"
esmat_files <- list.files(esmat_root, recursive = TRUE) %>%
  grep("ssGSEA_score.*\\.qs", ., value = TRUE)

type_pheno <- "survival"

bulks <- gsub(".*(GSE\\d+|TCGA_[A-Z]{4}).*", "\\1", esmat_files) %>% unique()

# ? read esmat from ssGSEA score files
purrr::walk(
  esmat_files,
  function(esmat_file_i) {
    # esmat score files under this bulk
    esmat_i <- qs::qread(file.path(esmat_root, esmat_file_i), nthreads = 2L)

    bulk <- gsub(".*(GSE\\d+|TCGA_[A-Z]{4}).*", "\\1", esmat_file_i) %>%
      unique()

    assign(paste0("esmat_", bulk), esmat_i, envir = .GlobalEnv)
  },
  .progress = TRUE
)

# # ? read test results for significance annotation
# diff_path <- "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/diff_test/survival/ov"

# test_files <- list.files(diff_path, recursive = TRUE) %>%
#   grep(
#     "csv$",
#     .,
#     value = TRUE
#   )

# purrr::walk(
#   test_files,
#   function(test_file_i) {
#     bulk <- gsub(".*(GSE\\d+|TCGA_[A-Z]{4}).*", "\\1", test_file_i)

#     assign(
#       paste0("test_res_", bulk),
#       data.table::fread(
#         file.path(diff_path, test_file_i)
#       ),
#       envir = .GlobalEnv
#     )
#   }
# )

# ? convert result to long format for visualization
purrr::walk(
  bulks,
  function(bulk_i) {
    # esmat score file paths under this bulk
    esmat_of_bulk_i <- qs::qread(
      file.path(esmat_root, esmat_files[grepl(bulk_i, esmat_files)]),
      nthreads = 2L
    )

    screen_method <- grepv(
      "(scissor|scPAS|scAB|scPP|DEGAS|LP_SGL|PIPET)$",
      colnames(esmat_of_bulk_i)
    ) %>%
      unique()

    esmat_of_bulk_i <- esmat_of_bulk_i %>%
      dplyr::select(
        dplyr::contains("ssGSEA"),
        dplyr::all_of(screen_method),
      ) %>%
      tidyr::pivot_longer(
        cols = tidyr::all_of(screen_method),
        names_to = "screen_method",
        values_to = "group"
      ) %>%
      tidyr::unite("cluster", screen_method, group) %>%
      dplyr::rename(
        pos_marker_score = dplyr::contains("pos"),
        neg_marker_score = dplyr::contains("neg")
      ) %>%
      tidyr::pivot_longer(
        cols = c(pos_marker_score, neg_marker_score),
        names_to = "type",
        values_to = "ssgsea_score"
      ) %>%
      dplyr::mutate(
        direction = dplyr::case_match(
          type,
          "pos_marker_score" ~ "Positive (↑)",
          "neg_marker_score" ~ "Negative (↓)"
        )
      )

    assign(
      paste0("esmat_", bulk_i),
      esmat_of_bulk_i,
      envir = .GlobalEnv
    )
  },
  .progress = TRUE
)

# ? plot boxplot for each bulk
library(ggplot2)
# library(ggsignif)
library(patchwork)
library(gghalves)

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
  rep(three_group_colors, 5)
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
  "LP_SGL_Neutral"
)


# # ? Plan A
# purrr::walk(
#   bulks,
#   function(bulk_i) {
#     esmat <- get(paste0("esmat_", bulk_i))
#     test_res_i <- get(paste0("test_res_", bulk_i))

#     p <- ggplot(
#       esmat,
#       aes(
#         x = cluster,
#         y = ssgsea_score,
#         fill = cluster, # box color
#       )
#     ) +
#       geom_boxplot(
#         outlier.alpha = 0.2,
#         outlier.size = 0.5,
#         outlier.colour = "#cececeff",
#         width = 0.65,
#         show.legend = FALSE
#       ) +
#       scale_fill_manual(
#         values = pallete
#       ) +
#       #   # 添加组间显著性连线与标签（关键！）
#       #   geom_signif(
#       #     comparisons = list(c("scAB_Other", "scAB_Positive")), # 指定比较的组名
#       #     map_signif_level = TRUE, # 自动映射 *** / ** / *
#       #     textsize = 6,
#       #   ) +
#       labs(
#         title = paste("ssGSEA Score of", bulk_i),
#         x = "Cluster",
#         y = "ssGSEA Score"
#       ) +
#       theme_minimal(base_size = 11.5) +
#       theme(
#         axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#         panel.grid.minor = element_blank()
#       )

#     obj_name <- paste0("plot_", gsub("[-.]", "_", bulk_i))
#     assign(obj_name, p, envir = .GlobalEnv)

#     cli::cli_alert_success("Created plot object: {obj_name}")
#   },
#   .progress = TRUE
# )

# ? Plan B
purrr::walk(
  bulks,
  function(bulk_i) {
    esmat <- get(paste0("esmat_", bulk_i))

    base_plot <- function(df_subset) {
      p <- ggplot(
        df_subset,
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
        scale_fill_manual(values = pallete) +
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
        )
      label_position <- df_subset %>%
        group_by(cluster) %>%
        summarise(
          y_pos = max(ssgsea_score, na.rm = TRUE),
          label = paste0("n=", n()),
          .groups = "drop"
        )

      p +
        ggrepel::geom_text_repel(
          data = label_position,
          aes(x = cluster, y = y_pos * 0.98 + 0.06, label = label),
          size = 3,
          fontface = "bold",
          box.padding = 0.5, # 文字周围的填充空间
          point.padding = 0.5, # 与数据点的距离
          segment.color = NA, # 不显示连接线（如需显示可去掉这行）
          direction = "y", # 主要沿 Y 轴方向排斥
          max.overlaps = Inf # 允许尝试所有重叠情况
        )
    }

    # Step 3: 拆分数据并绘图
    p_pos <- base_plot(
      esmat %>% dplyr::filter(direction == "Positive (↑)")
    ) +
      facet_grid("Positive Marker" ~ ., scales = "free_y", space = "free_y")

    p_neg <- base_plot(
      esmat %>% dplyr::filter(direction == "Negative (↓)")
    ) +
      facet_grid("Negative Marker" ~ ., scales = "free_y", space = "free_y")

    # Step 4: 拼接（重点：对齐 x 轴 & 共享图例/颜色）
    p_combined <- p_pos +
      p_neg +
      plot_layout(nrow = 1) +
      plot_annotation(
        title = paste("ssGSEA Score of", bulk_i),
        subtitle = paste0(
          "Bulk ",
          bulk_i,
          "; sc GSE165897; binary phenotype ",
          if (bulk_i == "GSE140082") {
            "High-grade/Low-grade"
          } else {
            "Malignant/LMP"
          }
        )
      )

    obj_name <- paste0("plot_", gsub("[-.]", "_", bulk_i))
    assign(obj_name, p_combined, envir = .GlobalEnv)

    cli::cli_alert_success("Created plot object: {obj_name}")
  },
  .progress = TRUE
)


# ? save boxplots
box_dir <- "boxplots"
if (!dir.exists(box_dir)) {
  dir.create(box_dir, recursive = TRUE)
}

purrr::walk(bulks, function(bulk_i) {
  ggsave(
    filename = file.path(
      box_dir,
      paste0("boxplot_ssGSEA_score_", bulk_i, ".png")
    ),
    plot = get(paste0("plot_", gsub("[-.]", "_", bulk_i))),
    width = 14,
    height = 7,
    dpi = 400
  )
  cli::cli_alert_success("Saved boxplot for {bulk_i}")
})
