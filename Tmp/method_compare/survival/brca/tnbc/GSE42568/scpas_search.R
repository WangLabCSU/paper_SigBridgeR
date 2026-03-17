library(dplyr)
library(data.table)
library(BiocParallel)

# ! BRCA
# ! sc- GSE161529
# ! bulk- GSE42568
# ! survival

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/survival/brca/tnbc/GSE42568"
)

data_path = "/home/data/sigbridger/benchmark_data/brca"
devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "seurat_tnbc.qs"))

bulk = qs::qread(
  file.path(data_path, "brca_bulkdata_GSE42568.qs")
)

pheno = qs::qread(file.path(
  data_path,
  "brca_pheno_GSE42568.qs"
))

surv_data = pheno %>%
  select("overall survival time_days:ch1", "overall survival event:ch1") %>%
  filter(`overall survival event:ch1` != "NA") %>% # cannot be changed to `!is.na()`
  rename("time" := 1, "status" := 2) %>%
  mutate_all(~ as.numeric(.))

bulk = bulk[, rownames(surv_data)]

if (!all(colnames(bulk) == rownames(surv_data))) {
  stop("bulk and surv_data not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}

# ! -------------------- scPAS ---------------------

# * random search, 50 times
set.seed(42568)
arg_samples <- data.frame(
  "nfeature" = sample(seq(500, 5000, by = 100), 50, replace = TRUE), # 第1维
  "imputation" = sample(c("None", "KNN", "ALRA"), 50, replace = TRUE), # 第2维
  "independent" = sample(c(FALSE, TRUE), 50, replace = TRUE) # 第3维
) %>%
  add_row(nfeature = 2000, imputation = "None", independent = TRUE)

save_path = 'scpas'
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# data.table::fwrite(
#     arg_samples,
#     file = file.path(save_path, "scpas_arg_samples.csv")
# )

split_arg_sample = arg_samples %>%
  mutate(group = (row_number() - 1) %/% 10 + 1) %>%
  group_split(group, .keep = FALSE)


# ! 并行设置
# mirai::daemons(2L)
# param <- MulticoreParam(workers = 2L, progressbar = TRUE)
# register(param)

setFuncOption(verbose = TRUE)

# * run scpas
for (j in c(2)) {
  # j = 2
  sub_arg_sample <- split_arg_sample[[j]]

  res_list <- lapply(
    seq_len(nrow(sub_arg_sample)),
    function(i) {
      nfeature_i = sub_arg_sample[i, "nfeature"][[1]]
      imputation_i = sub_arg_sample[i, "imputation"][[1]]
      independent_i = sub_arg_sample[i, "independent"][[1]]

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

      # 返回包含索引和结果的数据框
      return(data)
    }
    # ,
    # bulk = bulk,
    # seurat = seurat,
    # surv_data = surv_data,
    # sub_arg_sample = sub_arg_sample,
    # Screen = Screen,
    # glue = glue::glue
  )

  # *visualize
  gc()
  all_results <- do.call(cbind, res_list)
  rownames(all_results) = colnames(seurat)

  data.table::fwrite(
    all_results,
    file = file.path(save_path, paste0(j, "_scpas_random_search.csv")),
    row.names = TRUE
  )

  cli::cli_alert_success("Finished {j}th job")
}

# # ! ----------------可视化----------------

# # * cnv_status标记为肿瘤细胞，作为判断依据
# scpas_random_search = data.table::fread(
#     "scpas_random_search.csv",
# )

# seurat_tum = readRDS(
#     '/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_TNBCTum.rds'
# )
# tumor_cells = rownames(seurat_tum[[]])
# scpas_random_search$benchmark = ifelse(
#     colnames(seurat) %chin% tumor_cells,
#     TRUE,
#     FALSE
# )

# # * 计算指标
# source(
#     "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/compute_metrics.R"
# )
# metrics = compute_metrics(scpas_random_search)
# tmetrics = data.table::transpose(metrics)
# colnames(tmetrics) = rownames(metrics)
# rownames(tmetrics) = grep('process', colnames(scpas_random_search), value = T)

# arg_samples = as.data.table(arg_samples)
# arg_samples[, names(tmetrics) := tmetrics]

# data.table::fwrite(
#     arg_samples,
#     file = "scpas_arg_samples.csv",
#     row.names = TRUE
# )

# # * 图
# arg_samples <- arg_samples %>%
#     tidyr::unite("y_axis", imputation, independent, sep = " | ")

# library(ggplot2)
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
#     dpi = 400
# )

# # PID = 2861869
