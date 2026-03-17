# ! Negative Comparison
# ! Survival
# ! OV
# ! GSE140082
# ! GSE165897
library(dplyr)

setwd(
  '/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_negative_compare/survival/ov/GSE140082'
)

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_data/ov/GSE165897/GSE140082_ov_merged_seurat.qs",
  nthreads = 4L
)

score_path <- '/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_negative_compare/ov'
score_name <- paste0('ov_Sample_', 1:10, '_ssgsea_score')
score_files <- paste0(score_name, '.qs')

# ? read files
scores <- purrr::map2(
  score_files,
  score_name,
  function(file, name) {
    score_per_sample <- qs::qread(file.path(score_path, file))
    as.vector(score_per_sample)
  },
  .progress = TRUE
)
names(scores) <- score_name
scores_mat <- dplyr::bind_cols(scores)

meta <- seurat[[]]
screen_methods <- grepv(
  "sc[a-zA-Z]+$|LP_SGL$|PIPET$|DEGAS$",
  colnames(meta)
)
cli::cli_alert_info("screen_methods: {.field {screen_methods}}")
scores_mat <- cbind(scores_mat, meta[, screen_methods])

future::plan(future::multicore, workers = 5L)

run_ad_tests <- function(scores_mat, screen_methods, score_indices = 1:10) {
  # 确保输入为 data.frame 以便处理列名
  df <- as.data.frame(scores_mat)

  # 获取评分列名
  score_cols <- names(df)[score_indices]

  # 检查列名是否存在
  if (!all(screen_methods %in% names(df))) {
    stop("screen_methods 中的部分列名在 scores_mat 中不存在")
  }

  # 主循环：遍历每个分组方法
  purrr::map_dfr(screen_methods, function(g_col) {
    cli::cli_alert_info("Handling {.val {g_col}}")
    # 子循环：遍历每个评分指标
    furrr::future_map_dfr(
      score_cols,
      function(s_col) {
        # 构建临时数据框，去除 NA
        temp_df <- df %>%
          dplyr::select(dplyr::all_of(c(s_col, g_col))) %>%
          tidyr::drop_na()

        # 如果分组少于 2 个或样本太少，跳过
        if (length(unique(temp_df[[g_col]])) < 2 || nrow(temp_df) < 3) {
          cli::cli_abort("{s_col} {g_col} 分组少于 2 个或样本太少")
        }

        # 执行 AD 检验
        test_res <- tryCatch(
          {
            # kSamples::ad.test 支持公式接口
            kSamples::ad.test(
              reformulate(g_col, response = s_col),
              data = temp_df
            )
          },
          error = function(e) {
            cli::cli_warn("{e$message}")
            return(NULL)
          }
        )

        # 解析结果
        if (is.null(test_res)) {
          return(tibble::tibble(
            score_name = s_col,
            group_method = g_col,
            n_samples = test_res$N,
            n_groups = test_res$k,
            AS_stat = NA_real_,
            T_AV_stat = NA_real_,
            p_value = NA_real_,
            status = "Test Failed",
            test_name = test_res$test.name,
            n_ties = test_res$n.ties,
            sig = test_res$sig,
            warning = test_res$warning,
            null_dist1 = test_res$null.dist1,
            null_dist2 = test_res$null.dist2,
            method = test_res$method,
            n_sim = test_res$Nsim
          ))
        }
        tibble::tibble(
          score_name = s_col,
          group_method = g_col,
          n_samples = test_res$N,
          n_groups = test_res$k,
          AS_stat = test_res$ad[2, 1], # version 2
          T_AV_stat = test_res$ad[2, 2],
          p_value = test_res$ad[2, 3],
          status = "Success",
          test_name = test_res$test.name,
          n_ties = test_res$n.ties,
          sig = test_res$sig,
          warning = test_res$warning,
          null_dist1 = test_res$null.dist1,
          null_dist2 = test_res$null.dist2,
          method = test_res$method,
          n_sim = test_res$Nsim
        )
      },
      .progress = TRUE
    )
  })
}

ad_results <- run_ad_tests(
  scores_mat = scores_mat,
  screen_methods = screen_methods,
  score_indices = 1:10
)

data.table::fwrite(ad_results, 'ad_results.csv')
# p值 > 显著性水平（如0.05），则不能拒绝分布相同的原假设。
