library(dplyr)
library(data.table)
library(furrr)

# ! BRCA
# ! sc- GSE161529
# ! bulk- GSE42568
# ! survival

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/survival/brca/tnbc/GSE162228"
)

data_path = "/home/data/sigbridger/benchmark_data/brca"
devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "seurat_tnbc.qs"))

bulk = qs::qread(
  file.path(data_path, "brca_bulkdata_GSE162228.qs")
)

pheno = qs::qread(file.path(
  data_path,
  "brca_pheno_GSE162228.qs"
))

surv_data = pheno %>%
  select("overall survival (years):ch1", "characteristics_ch1.5") %>%
  rename("time" := 1) %>%
  mutate(
    status = case_when(
      characteristics_ch1.5 == "alive: Death" ~ 1,
      characteristics_ch1.5 == "alive: Alive" ~ 0
    )
  ) %>%
  select(-"characteristics_ch1.5") %>%
  mutate_all(~ as.numeric(.))

# bulk = bulk[, rownames(surv_data)]

if (!all(colnames(bulk) == rownames(surv_data))) {
  stop("bulk and surv_data not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}

# ! -------------------- scPAS ---------------------

# * random search, 100 times
set.seed(162228)
arg_samples <- data.frame(
  "nfeature" = sample(seq(500, 5000, by = 100), 50, replace = TRUE), # 第1维
  "imputation" = sample(c("None", "KNN", "ALRA"), 50, replace = TRUE), # 第2维
  "independent" = sample(c(FALSE, TRUE), 50, replace = TRUE) # 第3维
) %>%
  add_row(nfeature = 2000, imputation = "None", independent = TRUE)

# ! 并行设置
future::plan('multisession', workers = 2L)
options(future.globals.maxSize = 1024^3 * 4)

# * run scpas
if (.Platform$OS.type == "unix") {
  res_list <- furrr::future_map(
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
    },
    .options = furrr::furrr_options(
      seed = 123,
      packages = c("SigBridgeR", "furrr", "dplyr", "glue", 'future'),
      globals = c("seurat", "bulk", "surv_data", "arg_samples")
    ),
    .progress = TRUE
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

# 2004793

# ! ----------------可视化----------------

# * cnv_status标记为肿瘤细胞，作为判断依据
scpas_random_search = data.table::fread(
  "scpas_random_search.csv",
)

seurat_tum = readRDS(
  '/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_TNBCTum.rds'
)
tumor_cells = rownames(seurat_tum[[]])
scpas_random_search$benchmark = ifelse(
  colnames(seurat) %chin% tumor_cells,
  TRUE,
  FALSE
)

# * 计算指标
source(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/compute_metrics.R"
)
metrics = compute_metrics(scpas_random_search)
tmetrics = data.table::transpose(metrics)
colnames(tmetrics) = rownames(metrics)
rownames(tmetrics) = grep('process', colnames(scpas_random_search), value = T)

arg_samples = as.data.table(arg_samples)
arg_samples[, names(tmetrics) := tmetrics]

data.table::fwrite(
  arg_samples,
  file = "scpas_arg_samples.csv",
  row.names = TRUE
)

# * 图
arg_samples <- arg_samples %>%
  tidyr::unite("y_axis", imputation, independent, sep = " | ")

library(ggplot2)
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
  dpi = 400
)
