setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(data.table)

esmat_root <- "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/normalized_esmat/survival/ov"
esmat_files <- list.files(esmat_root, recursive = TRUE) %>%
  grep("ssGSEA_score_z.qs", ., value = TRUE)

bulks <- gsub(".*(GSE\\d+|TCGA_[A-Z]{4}).*", "\\1", esmat_files) %>% unique()

type_pheno <- "survival"

# ? read esmat from ssGSEA score files
purrr::walk(
  bulks,
  function(bulk_i) {
    # esmat score files under this bulk
    esmat_of_bulk_i <- esmat_files[grepl(bulk_i, esmat_files)]

    # screen methods under this bulk data
    screen_method <- gsub(
      ".*(scissor|scPAS|scAB|scPP).*",
      "\\1",
      esmat_of_bulk_i
    ) %>%
      unique()

    read_res <- purrr::walk2(
      screen_method,
      esmat_of_bulk_i,
      function(method_i, file_path) {
        assign(
          paste0("esmat_", bulk_i, "_", method_i), # data.table
          qs::qread(
            file.path(
              esmat_root,
              file_path
            ),
            nthreads = 4L
          ),
          envir = .GlobalEnv
        )
      }
    )
  },
  .progress = TRUE
)

# ? read test results for significance annotation
diff_path <- "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/normalized_diff_test/survival/ov"

test_files <- list.files(diff_path, recursive = TRUE) %>%
  grep(
    "csv$",
    .,
    value = TRUE
  )

purrr::walk(
  test_files,
  function(test_file_i) {
    bulk <- gsub(".*(GSE\\d+|TCGA_[A-Z]{4}).*", "\\1", test_file_i)

    assign(
      paste0("test_res_", bulk),
      data.table::fread(
        file.path(diff_path, test_file_i)
      ),
      envir = .GlobalEnv
    )
  }
)

# ? combine results
purrr::walk(
  bulks,
  function(bulk_i) {
    # esmat score file paths under this bulk
    esmat_of_bulk_i <- esmat_files[grepl(bulk_i, esmat_files)]

    screen_method <- gsub(
      ".*(scissor|scPAS|scAB|scPP).*",
      "\\1",
      esmat_of_bulk_i
    ) %>%
      unique()

    # combined object
    esmat_combined <- tibble::tibble()

    purrr::walk(
      screen_method,
      function(method_i) {
        esmat <- get(paste0("esmat_", bulk_i, "_", method_i))

        esmat_combined <<- dplyr::bind_rows(esmat_combined, esmat)
      }
    )

    assign(
      paste0("esmat_", bulk_i),
      esmat_combined,
      envir = .GlobalEnv
    )
  },
  .progress = TRUE
)


# ? plot boxplot for each bulk
library(ggplot2)
# library(ggsignif)
library(gghalves)

group_colors <- c(
  "Positive" = "#ff3333",
  "Negative" = "#386c9b",
  "Neutral" = "#CECECE"
)
pallete <- c(
  "#ff3333",
  "#CECECE",
  rep(group_colors, 3)
)
names(pallete) <- c(
  "scAB_Positive",
  "scAB_Other",
  "scPAS_Positive",
  "scPAS_Negative",
  "scPAS_Neutral",
  "scPP_Positive",
  "scPP_Negative",
  "scPP_Neutral",
  "scissor_Positive",
  "scissor_Negative",
  "scissor_Neutral"
)

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

purrr::walk(
  bulks,
  function(bulk_i) {
    esmat <- get(paste0("esmat_", bulk_i))
    test_res_i <- get(paste0("test_res_", bulk_i))

    p <- ggplot(
      esmat,
      aes(
        x = cluster,
        y = z_ssgsea_score,
        fill = cluster, # box color
      )
    ) +
      geom_half_boxplot(
        outlier.alpha = 0.2,
        outlier.size = 0.5,
        outlier.colour = "#cececeff",
        width = 0.65,
        alpha = 0.8,
        side = "l", # "l" = left
        show.legend = FALSE
      ) +
      geom_half_violin(
        side = "r", # "r" = right
        trim = FALSE, # 是否裁剪至数据范围（推荐 TRUE）
        alpha = 0.5, # 可调半透明度，避免遮挡 box
        show.legend = FALSE,
        width = 0.65,
        scale = "width"
      ) +
      scale_fill_manual(
        values = pallete
      ) +
      scale_y_continuous(
        breaks = scales::breaks_width(0.2), # 主刻度：0.0, 0.2, 0.4, ...
        minor_breaks = scales::breaks_width(0.1) # 次刻度：0.1, 0.3, 0.5, ...
      ) +

      labs(
        title = paste("ssGSEA Score of", bulk_i),
        subtitle = paste0(
          "Bulk ",
          bulk_i,
          "; sc GSE165897; survival as phenotype"
        ),
        x = "Cluster",
        y = "ssGSEA Score"
      ) +
      cowplot::theme_cowplot() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        panel.grid.minor = element_blank()
      )

    obj_name <- paste0("plot_", gsub("[-.]", "_", bulk_i))
    assign(obj_name, p, envir = .GlobalEnv)

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
      paste0("boxplot_z_ssGSEA_score_", bulk_i, ".png")
    ),
    plot = get(paste0("plot_", gsub("[-.]", "_", bulk_i))),
    width = 6,
    height = 5,
    dpi = 400
  )
  cli::cli_alert_success("Saved boxplot for {bulk_i}")
})
