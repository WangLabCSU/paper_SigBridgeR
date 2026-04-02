library(dplyr)
library(data.table)
library(Seurat)
# library(ggplot2)
library(BiocParallel)
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
seurat = qs::qread(file.path(data_path, "hgsoc_GSE165897_seurat.qs"))

bulk = qs::qread(
  file.path(data_path, "ov_bulkdata_GSE140082.qs"),
  nthreads = 4
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


# ! -------------------- scPP ---------------------

# * random search, 100 times
set.seed(123)
arg_samples <- data.frame(
  Log2FC_cutoff = sample(seq(0.05, 0.9, by = 0.05), 50, replace = TRUE),
  probs = sample(seq(0.05, 0.9, by = 0.05), 50, replace = TRUE)
) %>%
  add_row(Log2FC_cutoff = 0.585, probs = 0.2)


# * run scpp
if (.Platform$OS.type == "unix") {
  res_list <- bplapply(
    seq_len(nrow(arg_samples)),
    function(i) {
      tryCatch(
        {
          Log2FC_cutoff_i = arg_samples[i, 1]
          probs_i = arg_samples[i, 2]

          scpp_result = Screen(
            matched_bulk = bulk,
            sc_data = seurat,
            phenotype = surv_data,
            label_type = glue::glue("OS (M)_survival_{i}"),
            phenotype_class = "survival",
            screen_method = "scPP",
            # ref_group = 0, # actually, it's not used in survival data
            Log2FC_cutoff = Log2FC_cutoff_i,
            # estimate_cutoff = 0.2, # actually, it's not used in survival data
            probs = probs_i
          )

          pos_cell = (scpp_result$scRNA_data$scPP == "Positive")

          data = data.frame(
            pos_cell = pos_cell
          )
          colnames(data) = glue::glue("process_{i}")
          gc()

          # 返回包含索引和结果的数据框
          return(data)
        },
        error = function(e) {
          cli::cli_alert_danger(e$message)
          data = data.frame(
            pos_cell = FALSE
          )
          colnames(data) = glue::glue("process_{i}")
          return(data)
        }
      )
    }
  )
}

# *visualize
gc()
all_results <- do.call(cbind, res_list)
rownames(all_results) = colnames(seurat)

data.table::fwrite(
  all_results,
  file = "scpp_random_search.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("scPP random search completed."))

# ! ----- 可视化 ----

# * EOC细胞标记为肿瘤细胞，作为判断依据，注意有肿瘤细胞转移的情况
scpp_random_search = data.table::fread(
  "scpp_random_search.csv",
)
scpp_random_search$benchmark = setNames(
  grepl("EOC", seurat$cell_subtype),
  colnames(seurat)
)

# * 计算指标
source(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/ComputeMetrics.R"
)
metrics = ComputeMetrics(scpp_random_search)
t_metrics <- data.table::transpose(metrics)
colnames(t_metrics) <- rownames(metrics)
arg_samples = cbind(arg_samples, t_metrics)
data.table::fwrite(
  arg_samples,
  file = "scpp_arg_samples.csv",
  row.names = TRUE
)

library(ggplot2)

# * 图
p <- ggplot(
  arg_samples,
  aes(
    x = probs,
    y = Log2FC_cutoff,
    fill = F1
  )
) +
  geom_point(size = 6, alpha = 0.9, shape = 21, color = "black") +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "F1"
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = seq(0, 1, by = 0.1)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = seq(0, 1, by = 0.1)
  ) +
  labs(
    title = "Validation of the Screening Efficiency of scPP under Random Parameters",
    subtitle = "x = probs, y = Log2FC_cutoff",
    x = "probs",
    y = "Log2FC_cutoff"
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
  filename = "scpp_acc.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)

# PID=2352547
