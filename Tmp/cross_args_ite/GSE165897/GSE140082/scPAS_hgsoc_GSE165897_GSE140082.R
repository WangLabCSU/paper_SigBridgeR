library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(BiocParallel)

# !HGSOC.
# ! sc- GSE165897
# ! bulk- GSE140082
# ! survival

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/GSE165897/GSE140082"
)

data_path = "/home/data/sigbridger/benchmark_data/ov"
devtools::document("~/R/Project/R_code/SigBridgeR")

# ! 并行设置
param <- MulticoreParam(workers = 2)
register(param) # 注册为默认后端


# * load data
seurat = qs::qread(
  file.path(data_path, "hgsoc_GSE165897_seurat.qs"),
  nthreads = 4
)

bulk = qs::qread(
  file.path(data_path, "ov_bulkdata_GSE140082.qs"),
  nthreads = 4
)
# GSM4153781 has NA values, find median and replace NA
median = median(bulk[, "GSM4153781"], na.rm = TRUE)
na_indices <- which(is.na(bulk[, "GSM4153781"]))
bulk[na_indices, "GSM4153781"] <- median

pheno = qs::qread(file.path(data_path, "ov_pheno_GSE140082.qs"), nthreads = 4)

surv_data = select(
  pheno,
  "final_ostm:ch1",
  "final_osid:ch1"
) %>%
  rename("time" = 1, "status" = 2) %>%
  mutate_all(~ as.numeric(.))

# check
all(colnames(bulk) == rownames(surv_data))


# ! -------------------- scPAS ---------------------

# * random search, 100 times
set.seed(123)
arg_samples <- data.frame(
  "nfeature" = sample(seq(500, 5000, by = 100), 50, replace = TRUE), # 第1维
  "imputation" = sample(c("None", "KNN", "ALRA"), 50, replace = TRUE), # 第2维
  "independent" = sample(c(FALSE, TRUE), 50, replace = TRUE) # 第3维
) %>%
  add_row(nfeature = 2000, imputation = "None", independent = TRUE)


# * run scpas
if (.Platform$OS.type == "unix") {
  res_list <- bplapply(
    seq_len(nrow(arg_samples)),
    function(i) {
      nfeature_i = arg_samples[i, "nfeature"]
      imputation_i = arg_samples[i, "imputation"]
      independent_i = arg_samples[i, "independent"]

      if (imputation_i == "None") {
        scpas_result = Screen(
          matched_bulk = bulk,
          sc_data = seurat,
          phenotype = surv_data,
          label_type = glue::glue("OS (M)_survival_{i}"),
          phenotype_class = "survival",
          screen_method = "scPAS",
          alpha = NULL, # self-search
          independent = independent_i,
          imputation = FALSE,
          nfeature = nfeature_i,
          # assay = 'RNA',
          # network_class = "SC",
          # permutation_times = 2000,
          # FDR.threshold = 0.05
        )
      } else {
        scpas_result = Screen(
          matched_bulk = bulk,
          sc_data = seurat,
          phenotype = surv_data,
          label_type = glue::glue("OS (M)_survival_{i}"),
          phenotype_class = "survival",
          screen_method = "scPAS",
          alpha = NULL, # self-search
          independent = independent_i,
          imputation = TRUE,
          imputation_method = imputation_i,
          nfeature = nfeature_i,
          # assay = 'RNA',
          # network_class = "SC",
          # permutation_times = 2000,
          # FDR.threshold = 0.05
        )
      }

      pos_cell = (scpas_result$scRNA_data$scPAS == "Positive")

      data = data.frame(
        pos_cell = pos_cell
      )
      colnames(data) = glue::glue("process_{i}")
      gc()

      # 返回包含索引和结果的数据框
      return(data)
    }
  )
}

# *visualize
gc()
all_results <- do.call(cbind, res_list)
rownames(all_results) = colnames(seurat)

data.table::fwrite(
  all_results,
  file = "scpas_random_search.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("scPAS random search completed."))


# ! ----可视化----

# * EOC细胞标记为肿瘤细胞，作为判断依据，注意有肿瘤细胞转移的情况
scpas_random_search = data.table::fread(
  "scpas_random_search.csv",
)
scpas_random_search[, 1] <- NULL
scpas_random_search$benchmark = setNames(
  grepl("EOC", seurat$cell_subtype),
  colnames(seurat)
)

# * 计算指标
source(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/ComputeMetrics.R"
)
metrics = ComputeMetrics(scpas_random_search)
t_metrics <- data.table::transpose(metrics)
colnames(t_metrics) <- rownames(metrics)
arg_samples = cbind(arg_samples, t_metrics)

data.table::fwrite(
  arg_samples,
  file = "scpas_arg_samples.csv",
  row.names = TRUE
)

# * 图
arg_samples <- arg_samples %>%
  tidyr::unite("y_axis", imputation, independent, sep = " | ")

p <- ggplot(arg_samples, aes(x = nfeature, y = y_axis, fill = F1)) +
  geom_point(size = 6, alpha = 0.9, shape = 21, color = "black") +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "F1"
  ) +
  # 1. 保证横轴按真实数值排序，并强制每 500 一个刻度
  scale_x_continuous(
    breaks = seq(500, 5000, by = 250),
    labels = seq(500, 5000, by = 250)
  ) +
  labs(
    title = "Validation of the Screening Efficiency of scPAS under Random Parameters",
    subtitle = "x = nfeature, y = imputation | independent",
    x = "Number of features",
    y = NULL
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
  filename = "scpas_acc.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)

# PID = 1428120
