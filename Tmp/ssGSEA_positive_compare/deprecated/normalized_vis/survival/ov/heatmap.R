setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(tidyr)
library(data.table)

# ? read test results for significance annotation
diff_path <- "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/diff_test/survival/ov"

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

esmat_root <- "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/esmat/survival/ov"
esmat_files <- list.files(esmat_root, recursive = TRUE) %>%
  grep("ssGSEA_score.qs", ., value = TRUE)


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

# ? combine esmat according to bulk
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

# ? A function to calculate diff value
calculate_diff <- function(
  method_name = "scissor",
  group_name = "Positive",
  count = c(1, 3),
  diff = data.frame(),
  diff_grouped = data.frame()
) {
  this_cluster <- paste0(method_name, "_", group_name)

  all_cluster_of_method <- grepv(paste0("^", method_name, "_"), diff$cluster)

  rest_clusters <- setdiff(all_cluster_of_method, this_cluster)

  diff_value <- diff %>%
    filter(cluster == this_cluster) %>%
    pull(median_score)

  if (length(diff_value) == 0) {
    cli::cli_warn(
      "NA diff value for {this_cluster}, maybe this cluster does not exist. returning NA."
    )
    return(NA_real_)
  }

  if (count == 1) {
    # rest cluster is of length 1
    rest_value <- diff %>%
      filter(cluster == rest_clusters) %>%
      pull(median_score)
    return(diff_value[[1]] / rest_value[[1]])
  }

  #   count ==3
  possible_group <- c(
    paste0(rest_clusters[1], " + ", rest_clusters[2]),
    paste0(rest_clusters[2], " + ", rest_clusters[1])
  )

  rest_value <- diff_grouped %>%
    filter(
      method == method_name,
      pair_name %in% possible_group
    ) %>%
    pull(median_val)

  return(diff_value[[1]] / rest_value[[1]])
}


# ? calculate diff: median of group1 / median of rest
# ? add significance annotation
purrr::walk(
  bulks,
  function(bulk_i) {
    esmat_i <- get(paste0("esmat_", bulk_i))
    test_res_i <- get(paste0("test_res_", bulk_i))

    diff <- esmat_i %>%
      group_by(cluster) %>%
      summarise(median_score = median(ssgsea_score, na.rm = TRUE))

    diff_grouped <- esmat_i %>%
      mutate(
        method = sub("_.*", "", cluster)
      ) %>%
      group_by(method) %>%
      summarise(
        pairs = list(combn(unique(cluster), 2, simplify = FALSE)),
        .groups = "drop"
      ) %>%
      tidyr::unnest(pairs) %>%
      rowwise() %>%
      mutate(
        pair_name = paste(pairs[[1]], "+", pairs[[2]]),
        median_val = esmat_i %>%
          filter(cluster %in% pairs) %>%
          pull(ssgsea_score) %>% # ←←← 替换 `value` 为你的实际列名！
          median(na.rm = TRUE)
      ) %>%
      select(method, pair_name, median_val)

    diff_count <- diff_grouped %>%
      group_by(method) %>%
      summarise(
        count = n()
      )

    all_methods <- diff_count$method

    # scissor
    purrr::walk(
      all_methods,
      function(method_i) {
        purrr::walk(
          c("Positive", "Negative", "Neutral", "Other"),
          function(group_name) {
            if (
              method_i == "scAB" &&
                group_name %in% c("Negative", "Neutral") ||
                method_i != "scAB" && group_name == "Other"
            ) {
              return(NULL)
            }

            this_diff <- calculate_diff(
              method_name = method_i,
              group_name = group_name,
              count = diff_count$count[diff_count$method == method_i],
              diff = diff,
              diff_grouped = diff_grouped
            )

            test_res_i$diff[
              test_res_i$group == paste0(method_i, "_", group_name, " vs Rest")
            ] <<- this_diff
          }
        )
      },
      .progress = "Calculating diff values"
    )

    assign(
      paste0("test_res_", bulk_i),
      test_res_i %>%
        mutate(
          neg_log10_p = ifelse(p.value == 0, 300, -log10(p.value + 1e-300)), # 300 是上限
          group_simple = gsub(" vs Rest", "", group), # 去掉 " vs Rest"
          label = case_when(
            p.value < 1e-30 ~ "****",
            p.value < 1e-20 ~ "***",
            p.value < 1e-10 ~ "**",
            p.value < 0.01 ~ "*",
            TRUE ~ "NS"
          )
        ) %>%
        mutate(
          screen_group = paste0(screen_method, "\n", group_simple),
          screen_group = factor(screen_group, levels = unique(screen_group))
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

# # ? normalize diff for better visualization
combined <- combined %>%
  group_by(bulk) %>%
  mutate(
    normalized_diff = (diff - min(diff)) / (max(diff) - min(diff))
  )

# ? save combined results
data.table::fwrite(
  combined,
  file = "ov_test_data.csv"
)

cli::cli_alert_success(
  "Combined test results saved to {.file ov_test_data.csv}"
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
method_comparisons <- tribble(
  ~screen_method , ~comparison                ,
  "scAB"         , "scAB_Positive vs Rest"    ,
  "scAB"         , "scAB_Other vs Rest"       ,
  "scissor"      , "scissor_Positive vs Rest" ,
  "scissor"      , "scissor_Negative vs Rest" ,
  "scissor"      , "scissor_Neutral vs Rest"  ,
  "scPAS"        , "scPAS_Positive vs Rest"   ,
  "scPAS"        , "scPAS_Negative vs Rest"   ,
  "scPAS"        , "scPAS_Neutral vs Rest"    ,
  "scPP"         , "scPP_Positive vs Rest"    ,
  "scPP"         , "scPP_Negative vs Rest"    ,
  "scPP"         , "scPP_Neutral vs Rest"
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
  bulks,
  function(bulk_i) {
    # ? filter data for this bulk
    # ? because I found a wrong order when plotting all bulks together
    combined_i <- filter(combined, bulk == bulk_i)

    # ? complete missing combinations with NA
    plot_df <- select(
      combined_i,
      bulk,
      group,
      neg_log10_p,
      normalized_diff,
      label,
      screen_method
    ) %>%
      rename(
        comparison = group
      ) %>%
      right_join(method_comparisons, by = c("screen_method", "comparison")) %>%
      mutate(
        bulk = ifelse(is.na(bulk), bulk_i, bulk), # 补 bulk_i
        across(where(is.character), ~ tidyr::replace_na(.x, "NS")) # `label` is filled with "NS"
      )

    # 在调用 triangle_data 之前
    plot_df <- plot_df %>%
      mutate(
        # 将 screen_method 和 comparison 转换为因子，并指定你想要的顺序
        screen_method = factor(
          screen_method,
          levels = sort(c("scissor", "scPAS", "scPP", "scAB"))
        ),
        comparison = factor(
          comparison,
          levels = sort(plot_df$comparison)
        )
      ) %>%
      arrange(screen_method, comparison)

    plot_df2 <- triangle_data(plot_df, row = 1, col = 2)

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
        aes(upper.x, upper.y, fill = normalized_diff, group = group),
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
        limits = c(0, 1),
        na.value = "#e9e9e9ff", # ← NA 灰色
        name = "Norm.\nDiff."
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
        breaks = c(1:length(unique(plot_df2[[1]]))) + 0.5,
        expand = c(0, 0),
        labels = sort(unique(plot_df2[[1]]))
      ) +
      scale_y_continuous(
        expand = c(0, 0),
        breaks = c(1:length(unique(plot_df2[[2]]))) + 0.5,
        labels = sort(unique(plot_df2[[2]])),
        sec.axis = dup_axis()
      ) +
      labs(
        title = "OV ssGSEA Pos Ctrl",
        subtitle = glue::glue(
          "sc GSE165897 ;\nbulk {bulk_i};\nsurvival as phenotype"
        )
      ) +
      theme(
        axis.text.y.left = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(vjust = 0.5, size = 10, angle = 90),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10)
      )

    assign(paste0("p_", bulk_i), p, envir = .GlobalEnv)

    cli::cli_alert_success("Generated heatmap for bulk {.val {bulk_i}}")
  },
  .progress = TRUE
)


htmap_dir <- "heatmap"
if (!dir.exists(htmap_dir)) {
  dir.create(htmap_dir, recursive = TRUE)
}

# ? save plots
purrr::walk(bulks, function(bulk_i) {
  file_name <- glue::glue("ov_{bulk_i}_pos_ctrl_htmap.png")

  ggsave(
    filename = file.path(
      htmap_dir,
      file_name
    ),
    plot = get(paste0("p_", bulk_i)),
    width = 3.7,
    height = 6,
    dpi = 400
  )
  cli::cli_alert_success(
    "Heatmap saved to {.file {file_name}}"
  )
})
