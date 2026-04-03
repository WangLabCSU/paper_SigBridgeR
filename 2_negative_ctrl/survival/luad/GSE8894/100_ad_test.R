# ! Negative Comparison
# ! Survival
# ! luad
# ! GSE8894
# ! GSE123902
library(kSamples)
library(dplyr)
library(tidyr)

setwd(usethis::proj_path())

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/GSE8894/GSE8894_luad_merged_seurat.qs",
  nthreads = 4L
)

scores <- data.table::fread(
  "Tmp/ssGSEA_negative_compare/survival/luad/GSE8894/luad_sur_100reps_neg_ctrl_stat.csv"
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
          n_samples = NA_real_,
          n_groups = NA_real_,
          AS_stat = NA_real_,
          T_AV_stat = NA_real_,
          p_value = NA_real_,
          test_name = NA_character_,
          n_ties = NA_real_,
          sig = NA_real_,
          warning = NA_real_,
          null_dist1 = NA_real_,
          null_dist2 = NA_real_,
          method = NA_character_,
          n_sim = NA_real_,
          message = "分组数不足"
        ))
      }

      # 执行 k-sample Anderson-Darling 检验
      test_res <- ad.test(score_list)

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

data.table::fwrite(
  results,
  "Tmp/ssGSEA_negative_compare/survival/luad/GSE8894/rep100_ad_results.csv"
)
