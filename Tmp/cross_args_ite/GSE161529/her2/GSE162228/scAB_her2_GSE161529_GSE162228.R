library(dplyr)
library(Seurat)
library(BiocParallel)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/GSE161529/her2/GSE162228"
)

data_path = "/home/data/sigbridger/benchmark_data/brca"
devtools::document("~/R/Project/R_code/SigBridgeR")

# ! 并行设置
param <- MulticoreParam(workers = 2)
register(param) # 注册为默认后端

# * load data
seurat = qs::qread(file.path(data_path, "seurat_her2.qs"))

bulk = qs::qread(
  file.path(data_path, "brca_bulkdata_GSE162228.qs")
)
anyNA(bulk)

pheno = qs::qread(file.path(data_path, "brca_pheno_GSE162228.qs"))

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

# * random search, 100 times
set.seed(123)
alpha_samples <- data.frame(
  alpha = sample(seq(0.01, 1, 0.01), 50, replace = TRUE), # 第1维
  alpha2 = sample(seq(1e-08, 1e-03, 1e-08), 50, replace = TRUE) # 第2维
) %>%
  add_row(alpha = 0.005, alpha2 = 5e-05) # default parameters

# * run scab
# 确保输出目录存在
# if (!dir.exists("scAB_results")) {
#     dir.create("scAB_results")
# }

# * run scab with error handling
if (.Platform$OS.type == "unix") {
  res_list <- bplapply(
    seq_len(nrow(alpha_samples)),
    function(i) {
      tryCatch(
        {
          alpha1 <- alpha_samples[i, 1]
          alpha2 <- alpha_samples[i, 2]

          scab_result = Screen(
            matched_bulk = bulk,
            sc_data = seurat,
            phenotype = surv_data,
            label_type = glue::glue("OS (M)_survival_{i}"),
            phenotype_class = "survival",
            screen_method = "scAB",
            alpha = alpha1,
            alpha_2 = alpha2,
            # maxiter = 2000,
            # tred = 2
          )

          # qs::qsave(
          #     scab_result,
          #     file = glue::glue("scAB_results/scab_result_{i}.qs")
          # )

          pos_cell = (scab_result$scRNA_data$scAB == "Positive")

          data = data.frame(
            pos_cell = pos_cell
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
}

# 合并所有结果
all_results <- do.call(cbind, res_list)
rownames(all_results) = colnames(seurat)

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
alpha_samples$Accuracy = data.table::transpose(metrics)[[6]]

data.table::fwrite(
  alpha_samples,
  file = "scab_arg_samples.csv",
  row.names = TRUE
)
library(ggplot2)

plot = ggplot(alpha_samples, aes(x = alpha, y = alpha2, fill = Accuracy)) +
  geom_point(size = 6, alpha = .9, shape = 21, color = "black") +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "Accuracy"
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
