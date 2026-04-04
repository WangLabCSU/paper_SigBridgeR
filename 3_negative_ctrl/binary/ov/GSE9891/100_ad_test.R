# ! Negative Comparison
# ! Survival
# ! OV
# ! GSE9891
# ! GSE165897
library(kSamples)
library(dplyr)
library(tidyr)

setwd(usethis::proj_path())

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_binary/ov/GSE165897/GSE9891_hgsoc_merged_seurat.qs",
  nthreads = 4L
)

scores <- data.table::fread(
  "Tmp/ssGSEA_negative_compare/binary/ov/GSE9891/ov_bi_100reps_neg_ctrl_stat.csv"
)

scores <- tidyr::unite(
  scores,
  "ad_group",
  meta_column,
  label_value,
  sep = "_",
  remove = FALSE
)


# 按样本分组，对每个样本内的 ad_group 进行 ad.test
results <- scores %>%
  group_by(meta_column) %>%
  group_modify(
    ~ {
      # 将每个 ad_group 的 eamn_score 提取为列表
      score_list <- split(.x$mean_score, .x$ad_group)

      # 确保至少有 2 个分组
      if (length(score_list) < 2) {
        return(tibble(
          n_samples = NA,
          n_groups = NA,
          AS_stat = NA,
          T_AV_stat = NA,
          p_value = NA,
          test_name = NA,
          n_ties = NA,
          sig = NA,
          warning = NA,
          null_dist1 = NA,
          null_dist2 = NA,
          method = NA,
          n_sim = NA,
          message = "分组数不足"
        ))
      }

      # 执行 k-sample Anderson-Darling 检验
      test_res <- tryCatch(
        {
          ad.test(score_list)
        },
        error = function(e) {
          return(NULL)
        }
      )

      if (is.null(test_res)) {
        return(tibble(
          ad_statistic = NA,
          p_value = NA,
          n_groups = length(score_list),
          message = "检验失败"
        ))
      }

      tibble(
        n_samples = test_res$N,
        n_groups = test_res$k,
        AS_stat = test_res$ad[2, 1], # version 2
        T_AV_stat = test_res$ad[2, 2],
        p_value = test_res$ad[2, 3],
        test_name = test_res$test.name,
        n_ties = test_res$n.ties,
        sig = test_res$sig,
        warning = test_res$warning,
        null_dist1 = test_res$null.dist1,
        null_dist2 = test_res$null.dist2,
        method = test_res$method,
        n_sim = test_res$Nsim,
        message = "OK"
      )
    }
  ) %>%
  ungroup()

# 查看结果
print(results)

# 添加显著性标记
results <- results %>%
  mutate(
    significant = p_value < 0.05,
    sig_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

data.table::fwrite(results,"Tmp/ssGSEA_negative_compare/binary/ov/GSE9891/rep100_ad_results.csv")
