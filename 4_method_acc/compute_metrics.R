compute_metrics <- function(dt) {
  # library(data.table)
  # 获取列名
  all_cols <- names(dt)
  # 如果第一列是名称列，则计数大于2；T和F计数必定为2
  if (length(table(dt[, 1])) > 2) {
    process_cols <- all_cols[2:(ncol(dt) - 1)]
  } else {
    process_cols <- all_cols[1:(ncol(dt) - 1)]
  }

  # 提取benchmark列
  benchmark_logical <- as.logical(dt[['benchmark']])
  n_total <- length(benchmark_logical)

  # 计算基准统计
  pos_total <- sum(benchmark_logical, na.rm = TRUE)
  neg_total <- n_total - pos_total

  # 预计算基准向量
  actual_pos <- benchmark_logical
  actual_neg <- !benchmark_logical

  # 初始化结果向量
  n_proc <- length(process_cols)
  tpr_vec <- numeric(n_proc)
  fpr_vec <- numeric(n_proc)
  tnr_vec <- numeric(n_proc)
  fnr_vec <- numeric(n_proc)
  f1_vec <- numeric(n_proc)
  acc_vec <- numeric(n_proc)

  # 高效循环处理每个进程列
  for (i in seq_along(process_cols)) {
    pred_logical <- as.logical(dt[[process_cols[i]]])

    # 向量化计算混淆矩阵元素
    tp <- sum(pred_logical & actual_pos, na.rm = TRUE)
    fp <- sum(pred_logical & actual_neg, na.rm = TRUE)
    tn <- sum(!pred_logical & actual_neg, na.rm = TRUE)
    fn <- sum(!pred_logical & actual_pos, na.rm = TRUE)

    # 计算率值
    tpr_vec[i] <- ifelse(pos_total > 0, tp / pos_total, 0)
    fpr_vec[i] <- ifelse(neg_total > 0, fp / neg_total, 0)
    tnr_vec[i] <- ifelse(neg_total > 0, tn / neg_total, 0)
    fnr_vec[i] <- ifelse(pos_total > 0, fn / pos_total, 0)
    f1_vec[i] <- ifelse(tp + fp > 0, 2 * tp / (2 * tp + fp + fn), 0)
    acc_vec[i] <- (tp + tn) / (tp + fp + tn + fn)
  }

  # 构建结果data.table
  result <- data.table::data.table(
    TPR = tpr_vec,
    FPR = fpr_vec,
    TNR = tnr_vec,
    FNR = fnr_vec,
    F1 = f1_vec,
    Accuracy = acc_vec
  )
  result <- data.table::transpose(result)

  # 重命名以匹配进程列名
  data.table::setnames(result, seq_len(ncol(result)), process_cols)

  # 设置行名
  rownames(result) <- c("TPR", "FPR", "TNR", "FNR", "F1", "Accuracy")

  return(result)
}
