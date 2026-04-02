setwd(usethis::proj_path())

library(dplyr)

scores <- qs::qread(
  "Tmp/ssGSEA_negative_compare/brca/her2_Sample_100_ssgsea_score.qs",
  nthreads = 4L
)

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_binary/brca/HER2/GSE162228_brca_merged_seurat.qs",
  nthreads = 4L
)
meta <- seurat[[]]
screened_label <- meta[, grepv(
  "sc[a-zA-Z]+$|DEGAS$|LP_SGL$|PIPET$",
  colnames(meta)
)]

CalcGroupMeanScores <- function(
  score_mat,
  meta,
  meta_cols = NULL,
  sample_name = "Sample_1"
) {
  # 检查行名
  if (is.null(rownames(score_mat)) || is.null(rownames(meta))) {
    stop("score_mat 和 meta 必须包含行名 (cell IDs)。")
  }

  # 高效对齐：以 score_mat 的行顺序为基准，重排 meta
  # 使用 match 比 left_join 快得多
  common_cells <- intersect(rownames(score_mat), rownames(meta))
  if (length(common_cells) == 0) {
    stop("score_mat 和 meta 没有共同的细胞 ID。")
  }

  # 子集并对齐
  score_mat <- score_mat[common_cells, , drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]

  # 确定要处理的 meta 列
  if (is.null(meta_cols)) {
    meta_cols <- colnames(meta)
  } else {
    # 确保指定的列存在
    meta_cols <- intersect(meta_cols, colnames(meta))
    if (length(meta_cols) == 0) stop("指定的 meta_cols 在 meta 中不存在。")
  }

  # 2. 核心计算循环 (使用 rowsum 加速) ----------------------------------------
  results_list <- lapply(meta_cols, function(col_name) {
    group_vec <- meta[[col_name]]

    # 处理 NA 标签 (rowsum 会忽略 group 为 NA 的行，但为了统计计数需明确)
    # 这里我们暂时保留 NA 在 group_vec 中，rowsum 会自动排除它们，但 table 会统计它们
    # 为了结果整洁，我们显式移除 meta 中标签为 NA 的细胞参与计算
    valid_idx <- !is.na(group_vec)

    if (sum(valid_idx) == 0) {
      warning(paste("列", col_name, "全为 NA，跳过。"))
      return(NULL)
    }

    sub_mat <- score_mat[valid_idx, , drop = FALSE]
    sub_grp <- group_vec[valid_idx]

    # 确保 group 是字符型，防止因子水平问题
    sub_grp <- as.character(sub_grp)

    # --- 高效计算 ---
    # rowsum 计算每组每列的和 (速度极快)
    # ! 不是rowSum，是行聚合
    sum_mat <- rowsum(sub_mat, group = sub_grp, na.rm = TRUE)

    # table 计算每组计数
    count_vec <- as.vector(table(sub_grp))
    names(count_vec) <- rownames(sum_mat) # 确保名称对应

    # 计算均值
    # 注意：rowsum 返回的行名顺序通常与 table 一致，但为了安全最好匹配一下
    # 实际上 rowsum 和 table 对字符向量的排序逻辑一致 (sort unique)
    mean_mat <- sum_mat / count_vec

    # --- 格式化为长表格 ---
    # 此时数据量很小 (m2 * 标签数 * m)，转为 tibble 开销忽略不计
    res_df <- as.data.frame(mean_mat) %>%
      tibble::rownames_to_column(var = "label_value") %>%
      tidyr::pivot_longer(
        cols = -label_value,
        names_to = "score_group",
        values_to = "mean_score"
      ) %>%
      dplyr::mutate(
        meta_column = col_name,
        count = count_vec[match(label_value, names(count_vec))],
        log_count = log10(count + 1),
        sample = sample_name
      )

    return(res_df)
  })

  # 3. 合并结果 ---------------------------------------------------------------
  final_res <- dplyr::bind_rows(results_list)

  # 可选：整理列顺序，使其更接近您原函数的输出风格
  if (!is.null(final_res)) {
    final_res <- final_res %>%
      dplyr::select(
        meta_column,
        label_value,
        score_group,
        mean_score,
        count,
        log_count,
        sample
      )
  }

  return(final_res)
}

mean_scores <- CalcGroupMeanScores(
  scores,
  meta,
  meta_cols = colnames(screened_label),
  sample_name = "Sample_100_GSE162228"
)

dir_out <- "Tmp/ssGSEA_negative_compare/binary/brca/her2/GSE162228"

data.table::fwrite(
  mean_scores,
  file.path(dir_out, "her2_bi_100reps_neg_ctrl_stat.csv")
)

# mean_scores <- data.table::fread(file.path(
#   dir_out,
#   "luad_sur_100reps_neg_ctrl_stat.csv"
# ))
