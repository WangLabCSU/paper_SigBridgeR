library(dplyr)
library(Seurat)
library(BiocParallel)
# library(SigBridgeR)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/GSE165897/GSE140082"
)

devtools::document('/home/yyx/R/Project/R_code/SigBridgeR')
data_path = "/home/data/sigbridger/benchmark_data/ov"

# * load data
seurat = qs::qread(file.path(data_path, "hgsoc_GSE165897_seurat.qs"))

bulk = qs::qread(
  file.path(data_path, "ov_bulkdata_GSE140082.qs")
)
# GSM4153781 has NA values, find median and replace NA
median = median(bulk[, "GSM4153781"], na.rm = TRUE)
na_indices <- which(is.na(bulk[, "GSM4153781"]))
bulk[na_indices, "GSM4153781"] <- median

pheno = qs::qread(file.path(data_path, "ov_pheno_GSE140082.qs"))

surv_data = select(
  pheno,
  "final_ostm:ch1",
  "final_osid:ch1"
) %>%
  rename("time" = 1, "status" = 2) %>%
  mutate_all(~ as.numeric(.))

# check
all(colnames(bulk) == rownames(surv_data))


# ! -------------------- scAB ---------------------

rlogunif <- function(n, target, log_range = 1, digit = 4) {
  log_target <- log10(target) # median or center of distribution
  log_min <- log_target - log_range # 10^-log_range times smaller
  log_max <- log_target + log_range # 10^log_range times larger
  round(10^runif(n, log_min, log_max), digit)
}


# * random search, 100 times
set.seed(123)
arg_samples <- data.frame(
  alpha = rlogunif(50, 0.005), # 第1维
  alpha2 = rlogunif(50, 0.005), # 第2维
  maxiter = sample(seq(0, 10000, 100), 50, replace = TRUE),
  tred = sample(seq(0, 10, 1), 50, replace = TRUE)
) %>%
  add_row(alpha = 0.005, alpha2 = 0.005, maxiter = 2000, tred = 2) # default parameters


# * run scab with error handling
res_list <- bplapply(
  seq_len(nrow(arg_samples)),
  function(i) {
    rlang::try_fetch(
      {
        max_iter <- arg_samples$maxiter[i]
        tred <- arg_samples$tred[i]

        scab_result = Screen(
          matched_bulk = bulk,
          sc_data = seurat,
          phenotype = surv_data,
          label_type = glue::glue("OS (M)_survival_{i}"),
          phenotype_class = "survival",
          screen_method = "scAB",
          alpha = arg_samples$alpha, # select_alpha will be used
          alpha_2 = arg_samples$alpha2,
          maxiter = max_iter,
          tred = tred
        )

        # qs::qsave(
        #     scab_result,
        #     file = glue::glue("scAB_results/scab_result_{i}.qs")
        # )

        data = data.frame(
          pos_cell = (scab_result$scRNA_data$scAB == "Positive")
        )
        colnames(data) = glue::glue("process_{i}")

        # 返回包含索引和结果的数据框
        return(data)
      },
      error = function(e) {
        # 打印错误信息并返回包含索引和NA的数据框
        message(sprintf("Error at i=%d: %s", i, e$message))
        data = data.frame(pos_cell = FALSE)
        colnames(data) = glue::glue("process_{i}")
        return(data)
      }
    )
  }
)


# 合并所有结果
all_results <- do.call(cbind, res_list)
rownames(all_results) = colnames(seurat) # each cell is a row

data.table::fwrite(
  all_results,
  file = "scab_random_search.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("scAB random search completed."))

# # ! ----- Visualization ----

scab_random_search = data.table::fread(
  "scab_random_search.csv",
)
scab_random_search$benchmark = setNames(
  grepl("EOC", seurat$cell_subtype),
  colnames(seurat)
)

# * 计算指标
source(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/ComputeMetrics.R"
)
metrics = ComputeMetrics(scab_random_search)
t_metrics <- data.table::transpose(metrics)
colnames(t_metrics) <- rownames(metrics)
arg_samples = cbind(arg_samples, t_metrics)

data.table::fwrite(
  alpha_samples,
  file = "scab_arg_samples.csv",
  row.names = TRUE
)
library(ggplot2)

plot = ggplot(alpha_samples, aes(x = alpha, y = alpha2, fill = F1)) +
  geom_point(size = 6, alpha = 0.9, shape = 21, color = "black") +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "F1"
  ) +
  scale_color_gradient(low = "white", high = "red") +
  labs(
    title = "Validation of the Screening Efficiency of scAB under Random Parameters",
    subtitle = "x = alpha, y = alpha2",
    x = "alpha",
    y = "alpha2"
  ) +
  theme_minimal(base_size = 14) + # 全局字体基准
  theme(
    # 2. 轴文字放大
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    # 3. x 轴 45° 倾斜
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    # 4. 图例文字放大
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )

ggsave(
  filename = "scab_random_search.png",
  plot = plot,
  width = 10,
  height = 8,
  dpi = 300
)

# PID=3099497
