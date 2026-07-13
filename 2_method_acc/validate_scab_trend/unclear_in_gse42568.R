usethis::proj_activate(".")

her2_scab <- data.table::fread(
  "2_method_acc/brca_her2/arg_samples/scab_arg_samples2.csv" # ! another
)
tnbc_scab <- data.table::fread(
  "2_method_acc/brca_tnbc/arg_samples/scab_arg_samples2.csv" # !another
)

robust_kendall_perm <- function(
  x,
  y,
  alternative = "greater",
  B = 9999,
  seed = 123
) {
  set.seed(seed)
  n <- length(x)

  if (length(y) != n) {
    stop("x 和 y 的长度必须一致")
  }

  # 1. 计算观测到的 Kendall's tau_b (使用 pcaPP 的 C++ 底层加速，复杂度 O(N log N))
  obs_tau <- pcaPP::cor.fk(x, y)

  # 2. 置换检验核心循环
  # 预分配内存以提升 R 语言效率
  perm_taus <- numeric(B)

  for (i in 1:B) {
    # 打乱 y 的顺序，破坏 x 和 y 的配对关系，保留 y 的 Ties 结构
    y_perm <- sample(y)
    # 计算置换后的 tau
    perm_taus[i] <- pcaPP::cor.fk(x, y_perm)
  }

  # 3. 计算单侧/双侧 P 值
  # 使用 (k+1)/(B+1) 公式确保 p 值严格大于 0，避免 p=0 的尴尬
  if (alternative == "greater") {
    p_val <- (sum(perm_taus >= obs_tau) + 1) / (B + 1)
  } else if (alternative == "less") {
    p_val <- (sum(perm_taus <= obs_tau) + 1) / (B + 1)
  } else {
    p_val <- (sum(abs(perm_taus) >= abs(obs_tau)) + 1) / (B + 1)
  }

  # 4. 返回结果对象 (兼容 R 语言标准统计模型输出格式)
  res <- list(
    statistic = c(tau = obs_tau),
    p.value = p_val,
    method = paste("Kendall's tau_b with", B, "permutations"),
    alternative = alternative,
    null.value = c(tau = 0),
    perm_dist = perm_taus # 保留置换分布供后续可视化或诊断
  )
  class(res) <- "htest"
  return(res)
}


# 基于置换检验的 Kendall 相关性分析
her2_scab_tred_f1 <- robust_kendall_perm(
  x = her2_scab$tred,
  y = her2_scab$F1,
  seed = 123
)
her2_scab_tred_acc <- robust_kendall_perm(
  x = her2_scab$tred,
  y = her2_scab$Accuracy,
  seed = 124
)
# return to defaults
tnbc_scab_tred_f1 <- robust_kendall_perm(
  x = tnbc_scab$tred,
  y = tnbc_scab$F1,
  seed = 125
)
# return to defaults
tnbc_scab_tred_acc <- robust_kendall_perm(
  x = tnbc_scab$tred,
  y = tnbc_scab$Accuracy,
  seed = 126
)
all_scab_tred_f1 <- robust_kendall_perm(
  x = c(her2_scab$tred, tnbc_scab$tred),
  y = c(tnbc_scab$F1, her2_scab$F1),
  seed = 127
)
all_scab_tred_acc <- robust_kendall_perm(
  x = c(her2_scab$tred, tnbc_scab$tred),
  y = c(tnbc_scab$Accuracy, her2_scab$Accuracy),
  seed = 128
)
