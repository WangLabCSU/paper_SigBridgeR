pairwise_dist_test <- function(
  dt = data.table::data.table(),
  group_col = character(),
  value_col = character(),
  ...
) {
  # 获取所有方法名称
  groups <- unique(dt[, ..group_col]) %>% unlist()
  n_groups <- length(groups)

  # 初始化存储结果的列表
  results <- list()

  # 对所有配对进行循环
  cli::cli_progress_bar(
    name = "Test",
    type = "tasks",
    total = n_groups * (n_groups - 1) / 2
  )

  for (i in seq_len(n_groups - 1)) {
    for (j in (i + 1):n_groups) {
      g1 <- groups[i]
      g2 <- groups[j]

      # 提取两组数据（去除缺失值）
      x <- dt[dt[[group_col]] == g1, ..value_col] %>% unlist() %>% na.omit()
      y <- dt[dt[[group_col]] == g2, ..value_col] %>% unlist() %>% na.omit()

      # 跳过样本量太小的组（可选）
      if (length(x) < 3 || length(y) < 3) {
        next
      }

      # 执行 Anderson-Darling 两样本检验（k=2）
      # ad.test 要求输入为列表，每个元素是一个样本向量
      test_result <- tryCatch(
        expr = {
          l <- list(x, y)
          names(l) <- c(g1, g2)
          kSamples::ad.test(l)
        },
        error = function(e) NULL
      )

      # 记录结果
      results[[length(results) + 1]] <- if (!is.null(test_result)) {
        tibble::tibble(
          "test.name" = test_result$test.name,
          "k" = test_result$k,
          "ns.1" = test_result$ns[1],
          "ns.2" = test_result$ns[2],
          "N" = test_result$N,
          "n.ties" = test_result$n.ties,
          "sig" = test_result$sig,
          "ad.AD" = test_result$ad[2, 1],
          "ad.T.AD" = test_result$ad[2, 2],
          "AD.aympt. P-value" = test_result$ad[2, 3],
          "warning" = test_result$warning,
          "null.dist1" = test_result$null.dist1,
          "null.dist2" = test_result$null.dist2,
          "method" = test_result$method,
          "Nsim" = test_result$Nsim,
          group1 = g1,
          group2 = g2
        )
      } else {
        tibble::tibble(
          group1 = g1,
          group2 = g2,
          "warning" = TRUE
        )
      }

      cli::cli_progress_update()
    }
  }

  cli::cli_progress_done()
  dplyr::bind_rows(results)
}
