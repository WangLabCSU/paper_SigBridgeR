setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(tidyr)
library(data.table)

# ? read test results for significance annotation
diff_path <- "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/diff_test/binary/luad"

test_files <- list.files(diff_path, recursive = TRUE) %>%
  grep(
    "csv$",
    .,
    value = TRUE
  )

bulks <- gsub(".*(GSE\\d+|TCGA_[A-Z]{4}).*", "\\1", test_files) %>% unique()

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

esmat_root <- "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/esmat/binary/luad"
esmat_files <- list.files(esmat_root, recursive = TRUE) %>%
  grep("ssGSEA_score.*\\.qs", ., value = TRUE)


# ? read esmat from ssGSEA score files
purrr::walk(
  esmat_files,
  function(esmat_file_i) {
    bulk_i <- gsub(".*(GSE\\d+|TCGA_[A-Z]{4}).*", "\\1", esmat_file_i) %>%
      unique()

    assign(
      paste0("esmat_", bulk_i), # data.table
      qs::qread(
        file.path(
          esmat_root,
          esmat_file_i
        ),
        nthreads = 4L
      ),
      envir = .GlobalEnv
    )
  },
  .progress = TRUE
)

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

# ? calculate diff: mean of group1 / mean of rest
# ? add significance annotation
purrr::walk(
  bulks,
  function(bulk_i) {
    esmat_i <- get(paste0("esmat_", bulk_i))
    test_res_i <- get(paste0("test_res_", bulk_i))

    # for 0-1 scaling
    min_max <- esmat_i %>%
      dplyr::group_by(type) %>%
      dplyr::summarise(
        min_score = min(ssgsea_score),
        max_score = max(ssgsea_score)
      )
    esmat_i <- esmat_i %>%
      merge(min_max, by = "type") %>%
      mutate(
        scaled_ssgse_score = (ssgsea_score - min_score) /
          (max_score - min_score)
      )

    diff <- esmat_i %>%
      group_by(cluster, type) %>%
      summarise(mean_score = mean(scaled_ssgse_score, na.rm = TRUE))

    all_methods <- unique(gsub("_.*", "", diff$cluster))

    combined <- data.frame()

    # * rest median score
    purrr::walk(
      all_methods,
      function(method_i) {
        purrr::walk(
          c("Positive", "Negative", "Neutral", "Other"),
          function(group_name) {
            if (
              method_i %in%
                c("scAB", "DEGAS") &&
                group_name %in% c("Negative", "Neutral") ||
                !method_i %in% c("scAB", "DEGAS") &&
                  group_name == "Other"
            ) {
              return(NULL)
            }

            rest_mean <- esmat_i %>%
              dplyr::filter(
                grepl(method_i, cluster),
                !grepl(group_name, cluster) # exclude this group to get rest group
              ) %>%
              dplyr::group_by(type) %>%
              dplyr::summarise(rest_mean_score = mean(scaled_ssgse_score))

            mean_diff <- diff %>%
              dplyr::filter(cluster == paste0(method_i, "_", group_name))

            tb <- merge(mean_diff, rest_mean, by = "type")
            tb$diff <- tb$mean_score / tb$rest_mean_score

            combined <<- dplyr::bind_rows(combined, tb)
          }
        )
      },
      .progress = "Calculating diff values"
    )

    combined <- combined %>%
      dplyr::mutate(
        group = paste(gsub(".*_", "", cluster), "vs Rest"),
        method = gsub("_.*", "", cluster),
        type_sub = gsub("_.*", "", type)
      ) %>%
      tidyr::unite(
        "full_group",
        type_sub,
        method,
        group,
        sep = " | ",
        remove = FALSE
      )

    test_res_i <- test_res_i %>%
      dplyr::mutate(
        type_sub = gsub("_.*", "", `ssGSEA type`)
      ) %>%
      tidyr::unite(
        "full_group",
        type_sub,
        screen_method,
        group,
        sep = " | ",
        remove = FALSE
      )

    res <- merge(combined, test_res_i, by = "full_group")

    assign(
      paste0("test_res_", bulk_i),
      res %>%
        mutate(
          neg_log10_p = -log10(p.value + 1e-300), # 300 是上限
          label = case_when(
            p.value < 1e-3 ~ "***",
            p.value < 1e-2 ~ "**",
            p.value < 0.05 ~ "*",
            TRUE ~ "NS"
          )
        ) %>%
        mutate(
          group.x = factor(group.x, levels = unique(group.x)),
          group.y = factor(group.y, levels = unique(group.y))
        ),
      envir = .GlobalEnv
    )
    gc(verbose = FALSE)
  },
  .progress = TRUE
)

# ? combine all test results
combined <- data.frame()
for (bulk_i in bulks) {
  combined <- dplyr::bind_rows(combined, get(paste0("test_res_", bulk_i)))
}

# ? save combined results
data.table::fwrite(
  combined,
  file = "binary_luad_test_data.csv"
)

cli::cli_alert_success(
  "Combined test results saved to {.file binary_luad_test_data.csv}"
)

# combined <- data.table::fread("ov_test_data.csv")

# ----------------------------------------------------------------------------------------------------------------------
# ? plot heatmap, plan A
# library(ComplexHeatmap)

# df_plot <- select(combined, neg_log10_p, group, bulk) %>%
#   tidyr::pivot_wider(
#     names_from = group,
#     values_from = neg_log10_p
#   ) %>%
#   tibble::column_to_rownames("bulk") %>%
#   as.matrix()

# col_anno_method <- HeatmapAnnotation(
#   method = combined$screen_method[order(colnames(df_plot))],
#   col = list(
#     method = setNames(
#       c("#91bba8ff", "#e0c6a3ff", "#b8c4e9ff", "#c49fe2ff"),
#       c("scAB", "scissor", "scPAS", "scPP")
#     )
#   ),
#   height = unit(1, "cm"),
#   gap = unit(0.6, "mm")
# )

# Heatmap(
#   df_plot,
#   name = "-log10(P value)",
#   cluster_rows = FALSE,
#   cluster_columns = FALSE,
#   na_col = "#e9e9e9ff",
#   #   column_order = sort(colnames(df_plot)),
#   column_title = "Wilcox test group",
#   row_title = "Bulk",
#   rect_gp = gpar(col = "white", lwd = 2), # white cell border
#   top_annotation = col_anno_method,
#   bottom_annotation = HeatmapAnnotation(
#     empty = anno_empty(border = FALSE), # 让热图扁一点
#     height = unit(4, "mm")
#   )
# )

# ----------------------------------------------------------------------------------------------------------------------
# ? plot heatmap with ggplot2, plan B
library(ggplot2)


# ? all method comparisons, used to complete missing combinations
method_comparisons <- expand.grid(
  bulk = bulks,
  comparison = c(
    "scAB_Positive vs Rest",
    "scAB_Other vs Rest",
    "DEGAS_Positive vs Rest",
    "DEGAS_Other vs Rest",
    "PIPET_Positive vs Rest",
    "PIPET_Negative vs Rest",
    "PIPET_Neutral vs Rest",
    "scissor_Positive vs Rest",
    "scissor_Negative vs Rest",
    "scissor_Neutral vs Rest",
    "scPAS_Positive vs Rest",
    "scPAS_Negative vs Rest",
    "scPAS_Neutral vs Rest",
    "scPP_Positive vs Rest",
    "scPP_Negative vs Rest",
    "scPP_Neutral vs Rest",
    "LP_SGL_Positive vs Rest",
    "LP_SGL_Negative vs Rest",
    "LP_SGL_Neutral vs Rest"
  )
)

# ? generate triangle coordinates
triangle <- function(pairs, type = "up") {
  # 默认的上三角坐标基
  x = c(0, 0, 1)
  y = c(0, 1, 1)
  # 下三角的坐标基
  if (type == "lower") {
    x = c(0, 1, 1)
    y = c(0, 0, 1)
  }
  # 生成三角矩阵
  mat = do.call(
    rbind,
    apply(pairs, 1, function(row) {
      a = row[1]
      b = row[2]
      data.frame(
        x = x + a,
        y = y + b,
        group = paste(a, b, sep = "-")
      )
    })
  )
  return(mat)
}

# ? generate triangle data for ggplot2
triangle_data <- function(data, row = 1, col = 2) {
  # 保留原始因子水平顺序
  row_levels <- unique(data[[row]])
  col_levels <- unique(data[[col]])

  rows <- length(row_levels)
  cols <- length(col_levels)
  pairs <- merge(1:rows, 1:cols)

  # 获取上三角坐标
  upper <- triangle(pairs)
  colnames(upper) <- c(paste0("upper.", colnames(upper)[1:2]), "group")

  # 获取下三角坐标
  lower <- triangle(pairs, type = "lower")[1:2]
  colnames(lower) <- paste0("lower.", colnames(lower))

  # 合并坐标
  upper_lower <- bind_cols(upper, lower)
  pairs <- unite(data = pairs, col = "group", x, y, sep = "-")

  # **关键修改**：创建映射，使用原始顺序而非转换后的数值
  data$row_numeric <- as.numeric(factor(data[[row]], levels = row_levels))
  data$col_numeric <- as.numeric(factor(data[[col]], levels = col_levels))
  data$group <- paste(data$row_numeric, data$col_numeric, sep = "-")

  # 根据分组信息将坐标连接到数据中
  data <- data %>%
    right_join(upper_lower, by = "group")

  return(data)
}


purrr::walk(
  c("pos", "neg"),
  function(ssgsea_type) {
    # ? filter data for this bulk
    # ? because I found a wrong order when plotting all bulks together
    combined_i <- filter(combined, grepl(ssgsea_type, type_sub.y))

    # ? complete missing combinations with NA
    plot_df <- combined_i %>%
      unite(
        "comparison",
        screen_method,
        group.y,
        remove = FALSE
      ) %>%
      right_join(method_comparisons, by = c("bulk", "comparison")) %>%
      mutate(label = ifelse(is.na(label), "NA", label))

    # 在调用 triangle_data 之前
    plot_df <- plot_df %>%
      mutate(
        # 将 bulk 和 comparison 转换为因子，并指定你想要的顺序
        comparison = factor(
          comparison,
          levels = unique(sort(plot_df$comparison))
        ),
        bulk = factor(bulk, levels = unique(sort(plot_df$bulk)))
      ) %>%
      arrange(bulk, comparison)

    plot_df2 <- triangle_data(plot_df, row = "bulk", col = "comparison")

    # ? significance label
    points <- plot_df2
    points$x <- rep(0, nrow(points))
    points$y <- rep(0, nrow(points))

    # ? significance label position
    centers <- plot_df2 %>%
      group_by(bulk, comparison, group) %>%
      summarise(
        x_center = mean(lower.x),
        y_center = mean(lower.y),
        label = first(label),
        .groups = "drop"
      ) %>%
      filter(label != "" & !is.na(label))

    p <- ggplot(plot_df2) +
      geom_polygon(
        aes(upper.x, upper.y, fill = diff, group = group),
        colour = "grey",
        linewidth = 0.1
      ) +
      # diff颜色
      scale_fill_gradientn(
        colors = colorRampPalette(c(
          "#FFFFFF",
          "#FFED99",
          "#8ecde0ff",
          "#465d9bff"
        ))(10),
        limits = c(0.4, 2.5),
        na.value = "#e9e9e9ff", # ← NA 灰色
        name = "Diff"
      ) +
      ggnewscale::new_scale("fill") +
      # 显著性颜色
      geom_polygon(
        aes(lower.x, lower.y, fill = neg_log10_p, group = group),
        colour = "white",
        linewidth = 0.1
      ) +
      scale_fill_gradient(
        low = "#fceeeeff",
        high = "#d65456ff",
        limits = c(0, 300),
        na.value = "#e9e9e9ff",
        name = "-log10 (P value)"
      ) +
      geom_text(
        data = centers,
        aes(x = x_center + 0.1, y = y_center - 0.12, label = label),
        size = 2.4,
        fontface = "bold"
      ) +
      scale_x_continuous(
        breaks = c(1:length(unique(plot_df2[["bulk"]]))) + 0.5,
        expand = c(0, 0),
        labels = sort(unique(plot_df2[["bulk"]]))
      ) +
      scale_y_continuous(
        expand = c(0, 0),
        breaks = c(1:length(unique(plot_df2[["comparison"]]))) + 0.5,
        labels = sort(unique(plot_df2[["comparison"]])),
        sec.axis = dup_axis()
      ) +
      labs(
        title = "LUAD ssGSEA Pos Ctrl",
        subtitle = glue::glue(
          "sc GSE123902 ; {ssgsea_type} marker \nbinary phenotype\nTCGA_LUAD: Tumor/Normal"
        )
      ) +
      theme(
        axis.text.y.left = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(vjust = 0.5, size = 10, angle = 90),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10)
      )

    assign(paste0("p_", ssgsea_type), p, envir = .GlobalEnv)

    cli::cli_alert_success("Generated heatmap for bulk {.val {ssgsea_type}}")
  },
  .progress = TRUE
)


htmap_dir <- "heatmap"
if (!dir.exists(htmap_dir)) {
  dir.create(htmap_dir, recursive = TRUE)
}

# ? save plots
purrr::walk(c("pos", "neg"), function(ssgsea_type) {
  file_name <- glue::glue("luad_{ssgsea_type}marker_htmap.png")

  ggsave(
    filename = file.path(
      htmap_dir,
      file_name
    ),
    plot = get(paste0("p_", ssgsea_type)),
    width = 3.6,
    height = 8,
    dpi = 400
  )
  cli::cli_alert_success(
    "Heatmap saved to {.file {file_name}}"
  )
})
