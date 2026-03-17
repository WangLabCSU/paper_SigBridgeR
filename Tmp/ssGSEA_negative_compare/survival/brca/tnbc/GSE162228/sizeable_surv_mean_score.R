# ! Negative Comparison
# ! Survival
# ! BRCA - TNBC
# ! GSE162228
# ! GSE161529

setwd(
  '/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_negative_compare/survival/brca/tnbc/GSE162228'
)

score_path <- '/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_negative_compare/brca'
score_name <- paste0('tnbc_Sample_', 1:10, '_ssgsea_score')
score_files <- paste0(score_name, '.qs')

# ? read files
purrr::walk2(
  score_files,
  score_name,
  function(file, name) {
    assign(name, qs::qread(file.path(score_path, file)), envir = .GlobalEnv)
  },
  .progress = TRUE
)

# ? A function to calculate the mean score of each group
MeanScore = function(seurat_obj, es.mat, sample = "GSE", i = NULL) {
  library(tidyr)

  meta <- seurat_obj[[]] %>%
    tibble::rownames_to_column("cell")

  # 处理表达矩阵并合并元数据
  es.df <- Matrix::t(es.mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::left_join(meta, by = "cell")

  score_row_name = if (is.null(i)) {
    'Neg_ctrl_Sample'
  } else {
    paste0('Neg_ctrl_Sample_', i)
  }

  calculate_counts <- function(df, method_col) {
    if (method_col %in% colnames(df)) {
      df %>%
        dplyr::group_by(!!sym(method_col)) %>%
        dplyr::summarise(
          count = n(), # 计算每组的计数
          dplyr::across(
            score_row_name,
            mean,
            na.rm = TRUE
          )
        ) %>%
        pivot_longer(
          cols = score_row_name,
          names_to = "cluster",
          values_to = "mean_score"
        ) %>%
        pivot_longer(
          cols = all_of(method_col),
          names_to = "screen_method",
          values_to = "screen_result"
        )
    } else {
      NULL
    }
  }

  # 计算每个方法的计数和均值
  mean_scores_scissor <- calculate_counts(es.df, "scissor")
  mean_scores_scpas <- calculate_counts(es.df, "scPAS")
  mean_scores_scab <- calculate_counts(es.df, "scAB")
  mean_scores_scpp <- calculate_counts(es.df, "scPP")
  mean_scores_degas <- calculate_counts(es.df, "DEGAS")
  mean_scores_lp_sgl <- calculate_counts(es.df, "LP_SGL")

  # 合并所有数据
  mean_scores <- rbind(
    mean_scores_scissor,
    mean_scores_scpas,
    mean_scores_scab,
    mean_scores_scpp,
    mean_scores_degas,
    mean_scores_lp_sgl
  ) %>%
    unite("screen", screen_method, screen_result, sep = "_") %>%
    dplyr::mutate("Sample" = sample)

  # 对计数进行对数转换（避免log(0)的问题）
  mean_scores <- mean_scores %>%
    mutate(
      log_count = log10(count + 1), # +1 避免log(0)
      bubble_size = sqrt(count) # 或者使用平方根转换，视觉效果更好
    )
}

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_data/brca/TNBC/GSE162228_tnbc_merged_seurat.qs",
  nthreads = 4L
)

# ? Calculate the mean score of each group
mean_scores <- purrr::imap(
  score_name,
  function(name, index) {
    MeanScore(
      seurat_obj = seurat,
      es.mat = get(name),
      sample = 'GSE162228',
      i = index
    )
  },
  .progress = TRUE
)

names(mean_scores) <- score_name
tbl_mean_scores <- dplyr::bind_rows(mean_scores)

# ? Convert the cluster name to factor to make the order of the cluster fixed
tbl_mean_scores$cluster <- factor(
  tbl_mean_scores$cluster,
  levels = unique(tbl_mean_scores$cluster)
)


bubble <- ggplot2::ggplot(
  tbl_mean_scores,
  ggplot2::aes(x = cluster, y = screen, fill = mean_score, size = log_count)
) +
  ggplot2::geom_point(alpha = 0.9, shape = 21, color = "black") +
  ggplot2::scale_size_continuous(
    name = "Cell Count (log10+1)",
    range = c(3, 5), # 调整气泡大小范围
    breaks = log10(c(500, 1000, 5000, 10000, 50000) + 1), # 对数断点
    labels = c('500', "1K", '5K', "10K", "50K") # 原始数值标签
  ) +
  ggplot2::scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    name = "Mean Score"
  ) +
  ggplot2::labs(
    title = "Mean ssGSEA Scores by Group",
    subtitle = "sc-GSE161529 TNBC, bulk-GSE162228",
    y = "Mean ssGSEA Score",
    x = "Group"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

ggplot2::ggsave(
  filename = 'Sizable_Neg_ctrl_tnbc_GSE162228_mean_score.png',
  bubble,
  dpi = 400,
  width = 6,
  height = 6
)
