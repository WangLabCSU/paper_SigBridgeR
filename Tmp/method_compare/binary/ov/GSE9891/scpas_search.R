library(dplyr)
# library(tidyr)
library(data.table)

# library(BiocParallel)

# ! HGSOC.
# ! sc- GSE165897
# ! bulk- GSE9891
# ! binary

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/binary/ov/GSE9891"
)

data_path = "/home/data/sigbridger/benchmark_data/ov"
devtools::document("~/R/Project/R_code/SigBridgeR")

# ! 并行设置
# param <- MulticoreParam(workers = 2)
# register(param) # 注册为默认后端

# * load data
seurat = qs::qread(file.path(data_path, "hgsoc_GSE165897_seurat.qs"))

bulk = qs::qread(
  file.path(data_path, "ov_bulkdata_GSE9891.qs"),
  nthreads = 4
)

pheno = qs::qread(file.path(data_path, "ov_pheno_GSE9891.qs"))

pheno_bi = setNames(
  case_when(
    pheno$characteristics_ch1.1 == "Type : LMP" ~ 0,
    pheno$characteristics_ch1.1 == "Type : Malignant" ~ 1
  ),
  pheno$geo_accession
)

# ! 匹配
bulk = bulk[, names(pheno_bi)]

if (any(colnames(bulk) != names(pheno_bi))) {
  stop("bulk and pheno_bi not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}
if (anyNA(pheno_bi)) {
  stop("pheno_bi has NA")
}

# ! -------------------- scPAS ---------------------

# * random search, 100 times
set.seed(9891)
arg_samples <- data.frame(
  "nfeature" = sample(seq(500, 5000, by = 100), 50, replace = TRUE), # 第1维
  "imputation" = sample(c("None", "KNN", "ALRA"), 50, replace = TRUE), # 第2维
  "independent" = sample(c(FALSE, TRUE), 50, replace = TRUE) # 第3维
) %>%
  add_row(nfeature = 2000, imputation = "None", independent = TRUE)


# * run scpas
if (.Platform$OS.type == "unix") {
  res_list <- lapply(
    seq_len(nrow(arg_samples)),
    function(i) {
      nfeature_i = arg_samples[i, "nfeature"][[1]]
      imputation_i = arg_samples[i, "imputation"][[1]]
      independent_i = arg_samples[i, "independent"][[1]]

      if (imputation_i == "None") {
        scpas_result = Screen(
          matched_bulk = bulk,
          sc_data = seurat,
          phenotype = pheno_bi,
          label_type = glue::glue("Type"),
          phenotype_class = "binary",
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
          phenotype = pheno_bi,
          label_type = glue::glue("Type"),
          phenotype_class = "binary",
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

library(ggplot2)

# * EOC细胞标记为肿瘤细胞，作为判断依据，注意有肿瘤细胞转移的情况
scpas_random_search = data.table::fread(
  "scpas_random_search.csv",
)
scpas_random_search$benchmark = setNames(
  grepl("EOC", seurat$cell_subtype),
  colnames(seurat)
)

# * 计算指标
source(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/compute_metrics.R"
)
metrics = compute_metrics(scpas_random_search)
arg_samples <- bind_cols(
  arg_samples,
  data.table::transpose(metrics)
)
colnames(arg_samples) <- c(
  "nfeature",
  "imputation",
  "independent",
  rownames(metrics)
)

data.table::fwrite(
  arg_samples,
  file = "scpas_arg_samples.csv",
  row.names = TRUE
)

# * 图
arg_samples <- arg_samples %>%
  tidyr::unite("y_axis", imputation, independent, sep = " | ")

p <- ggplot(arg_samples, aes(x = nfeature, y = y_axis, fill = Accuracy)) +
  geom_point(size = 6, alpha = .9, shape = 21, color = "black") +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "Accuracy"
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

# # PID = 1490637
