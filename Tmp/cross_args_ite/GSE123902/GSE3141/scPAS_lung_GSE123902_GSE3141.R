library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(BiocParallel)

# ! LUAD
# ! sc - GSE123902
# ! bulk - GSE3141
# ! survival

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/GSE123902/GSE3141"
)


data_path = "/home/data/sigbridger/benchmark_data/lung"
devtools::document("~/R/Project/R_code/SigBridgeR")

setFuncOption(verbose = F)


# ! 并行设置
param <- BiocParallel::MulticoreParam(
  workers = 2,
  progressbar = TRUE,
  log = TRUE, # 启用日志
  threshold = "INFO" # 信息级别
)
register(param) # 注册为默认后端


# * load data
seurat = qs::qread(
  file.path(data_path, "luad_GSE123902_seurat.qs"),
  nthreads = 4
)


bulk = qs::qread(
  file.path(data_path, "lung_bulkdata_GSE3141.qs")
)
!anyNA(bulk)

pheno = data.table::fread(
  "/home/data/data-resource/single-cell/Lung_Cancer/LUAD/lung.cancer.adeno.gse3141.hgu133plus2_entrezcdf.tsv"
)

surv_data = pheno %>%
  tibble::column_to_rownames("Array") %>%
  rename("time" = 2, "status" = 1) %>%
  .[, c("time", "status")]
# check
# bulk = bulk[, rownames(surv_data)]
all(colnames(bulk) == rownames(surv_data))


# ! -------------------- scPAS ---------------------

# * random search, 100 times
set.seed(3141)
arg_samples <- data.frame(
  "nfeature" = sample(seq(500, 5000, by = 100), 50, replace = TRUE), # 第1维
  "imputation" = sample(c("None", "KNN", "ALRA"), 50, replace = TRUE), # 第2维
  "independent" = sample(c(FALSE, TRUE), 50, replace = TRUE) # 第3维
) %>%
  add_row(nfeature = 2000, imputation = "None", independent = TRUE)


# * run scpas
# ! 并行设置
param <- BiocParallel::MulticoreParam(
  workers = 2,
  progressbar = TRUE,
  log = TRUE, # 启用日志
  threshold = "INFO" # 信息级别
)
register(param) # 注册为默认后端

if (.Platform$OS.type == "unix") {
  res_list <- BiocParallel::bplapply(
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
          label_type = glue::glue("OS (days)_survival_{i}"),
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
        scpas_result = rlang::try_fetch(
          Screen(
            matched_bulk = bulk,
            sc_data = seurat,
            phenotype = surv_data,
            label_type = glue::glue("OS (days)_survival_{i}"),
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
          ),
          error = list(
            scRNA_data = list(scPAS = rep("Neutral", ncol(seurat)))
          )
        )
      }

      pos_cell = (scpas_result$scRNA_data$scPAS == "Positive")

      data = data.frame(
        pos_cell = pos_cell
      )
      col_name = glue::glue("process_{i}")
      colnames(data) = col_name
      message(col_name, " finished.")
      gc()

      # 返回包含索引和结果的数据框
      return(data)
    }
  )
}

# *visualize
gc()
all_results <- do.call(cbind, res_list)
all_results$cell_id = colnames(seurat)

vroom::vroom_write(
  all_results,
  file = "scpas_random_search.csv",
  delim = ","
)

vroom::vroom_write(arg_samples, "scpas_arg_samples.csv", delim = ",")

cli::cli_alert_success(crayon::green("scPAS random search completed."))

# ! ----可视化----

# # # ! HER2 seurat没有可用指标，舍弃
# scpas_random_search = data.table::fread(
#     "scpas_random_search.csv",
# )
# scpas_random_search = scpas_random_search %>%
#     tibble::column_to_rownames("cell_id") %>%
#     tibble::rownames_to_column("cell_id")
# scpas_random_search$benchmark = setNames(
#     grepl("tumor", seurat$cnv_status),
#     colnames(seurat)
# )

# # * 计算指标
# source(
#     "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/ComputeMetrics.R"
# )

# metrics = ComputeMetrics(scpas_random_search)
# arg_samples$Accuracy = data.table::transpose(metrics)[[6]]

# data.table::fwrite(
#     arg_samples,
#     file = "scpas_arg_samples.csv",
#     row.names = TRUE
# )

# # * 图
# arg_samples <- arg_samples %>%
#     tidyr::unite("y_axis", imputation, independent, sep = " | ")

# p <- ggplot(arg_samples, aes(x = nfeature, y = y_axis, fill = Accuracy)) +
#     geom_point(size = 6, alpha = .9, shape = 21, color = "black") +
#     scale_fill_gradient(
#         low = "white",
#         high = "red",
#         name = "Accuracy"
#     ) +
#     # 1. 保证横轴按真实数值排序，并强制每 500 一个刻度
#     scale_x_continuous(
#         breaks = seq(500, 5000, by = 250),
#         labels = seq(500, 5000, by = 250)
#     ) +
#     labs(
#         title = "Validation of the Screening Efficiency of scPAS under Random Parameters",
#         subtitle = "x = nfeature, y = imputation | independent",
#         x = "Number of features",
#         y = NULL
#     ) +
#     theme_minimal(base_size = 14) + # 全局字体基准
#     theme(
#         # 2. 轴文字放大
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 13),
#         # 3. x 轴 45° 倾斜
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         # 4. 图例文字放大
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 13)
#     )

# ggsave(
#     filename = "scpas_acc.png",
#     plot = p,
#     width = 10,
#     height = 8,
#     dpi = 300
# )

# # PID = 3917335
