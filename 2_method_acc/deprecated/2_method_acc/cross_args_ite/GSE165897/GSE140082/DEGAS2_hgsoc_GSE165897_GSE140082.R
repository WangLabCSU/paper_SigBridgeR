library(dplyr)
library(data.table)
library(Seurat)
library(ggplot2)
# library(BiocParallel)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/GSE165897/GSE140082"
)

data_path = "/home/data/sigbridger/benchmark_data/ov"
devtools::document("~/R/Project/R_code/SigBridgeR")

# # ! 并行设置
# param <- MulticoreParam(workers = 2)
# register(param) # 注册为默认后端
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

pheno_bi = pheno[
  pheno$`newgrade:ch1` != "NA" &
    !is.na(pheno$`newgrade:ch1`),
]
pheno_bi = setNames(
  case_when(
    pheno_bi$`newgrade:ch1` == "high.grade" ~ 1,
    pheno_bi$`newgrade:ch1` == "low.grade" ~ 0
  ),
  pheno_bi$geo_accession
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


# ! -------------------- degas ---------------------

# * random search, 100 times
set.seed(123)
arg_samples <- data.frame(
  lamb1 = sample(2:10, 50, replace = TRUE),
  lamb2 = sample(2:10, 50, replace = TRUE),
  lamb3 = sample(2:10, 50, replace = TRUE)
) %>%
  add_row(lamb1 = 3, lamb2 = 3, lamb3 = 3)

SigBridgeR::setThreads(
  8L,
  tf_config = list(
    intra_op = 8L,
    inter_op = 8L
  )
)


# * run DEGAS
res_list <- lapply(
  seq_len(nrow(arg_samples)),
  function(i) {
    tryCatch(
      {
        lamb1 <- arg_samples[i, "lamb1"][[1]]
        lamb2 <- arg_samples[i, "lamb2"][[1]]
        lamb3 <- arg_samples[i, "lamb3"][[1]]

        result <- Screen(
          matched_bulk = bulk,
          sc_data = seurat,
          phenotype = surv_data,
          label_type = glue::glue("OS (M)_survival_{i}"),
          phenotype_class = "survival",
          screen_method = "DEGAS",
          degas_params = list(
            DEGAS.lambda1 = lamb1,
            DEGAS.lambda2 = lamb2,
            DEGAS.lambda3 = lamb3
          )
        )

        pos = (result$scRNA_data$DEGAS == "Positive")

        data = data.frame(
          pos_cell = pos
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

# data <- qs::qread("/home/data/sigbridger/merged_without_pipet.qs")

# *visualize
gc()
all_results <- do.call(cbind, res_list)
rownames(all_results) = colnames(seurat)

data.table::fwrite(
  all_results,
  file = "degas2_random_search.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("DEGAS random search completed."))

# # ! ----- 可视化 ----

# * EOC细胞标记为肿瘤细胞，作为判断依据，注意有肿瘤细胞转移的情况
degas_random_search = data.table::fread(
  "degas2_random_search.csv",
)
degas_random_search$benchmark = setNames(
  grepl("EOC", seurat$cell_subtype),
  colnames(seurat)
)

for (i in seq_len(ncol(degas_random_search))) {
  set(
    degas_random_search,
    j = i,
    value = as.logical(degas_random_search[[i]] == "Positive")
  )
}

# * 计算指标
source(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/ComputeMetrics.R"
)
metrics = ComputeMetrics(degas_random_search)
t_metrics <- data.table::transpose(metrics)
colnames(t_metrics) <- rownames(metrics)
arg_samples = cbind(arg_samples, t_metrics)

data.table::fwrite(
  arg_samples,
  file = "degas2_arg_samples.csv",
  row.names = TRUE
)

# * 图
p <- ggplot(
  arg_samples,
  aes(
    x = bag_depth,
    y = ff_depth,
    fill = F1,
    shape = arch
  )
) +
  ggbeeswarm::geom_quasirandom(
    size = 6,
    alpha = 0.9,
    color = "black",
    method = "quasirandom", # 或 "swarm", "quasirandom"
    groupOnX = TRUE # 按x轴分组躲避
  ) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "F1"
  ) +
  labs(
    title = "Validation of the Screening Efficiency of degas under Random Parameters",
    subtitle = "x = bag_depth, y = ff_depth",
    x = "bag_depth",
    y = "ff_depth"
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
  filename = "degas2_acc.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)
