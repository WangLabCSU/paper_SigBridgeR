library(dplyr)
library(Seurat)
library(BiocParallel)
library(SigBridgeR)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(file.path(
  usethis::proj_path(),
  "Tmp/method_compare/binary/ov/GSE140082"
))
SigBridgeR::setThreads(10L)

# devtools::document('/home/yyx/R/Project/R_code/SigBridgeR')
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


pheno = pheno[
  pheno$`newgrade:ch1` != "NA" &
    !is.na(pheno$`newgrade:ch1`),
]
pheno = setNames(
  case_when(
    pheno$`newgrade:ch1` == "high.grade" ~ 1,
    pheno$`newgrade:ch1` == "low.grade" ~ 0
  ),
  pheno$geo_accession
)

# ! 匹配
bulk = bulk[, names(pheno)]


# check
all(colnames(bulk) == names(pheno))

rlogunif <- function(n, target, log_range = 1, digit = 4) {
  log_target <- log10(target) # median or center of distribution
  log_min <- log_target - log_range # 10^-log_range times smaller
  log_max <- log_target + log_range # 10^log_range times larger
  round(10^runif(n, log_min, log_max), digit)
}

distance_choices <- c(
  "cosine",
  "pearson",
  "spearman",
  "kendall",
  "euclidean",
  "maximum"
)

# * random search, 50 times
set.seed(123)
arg_samples <- data.frame(
  distance = sample(distance_choices, 50, replace = TRUE), # 第1维
  nPerm = sample(seq(500, 5000, 100), 50, replace = TRUE),
  log2FC = sample(seq(0.5, 2, 0.01), 50, replace = TRUE)
) %>%
  add_row(distance = "cosine", nPerm = 1000L, log2FC = 1L) # default parameters


# * run scab with error handling
res_list <- lapply(
  seq_len(nrow(arg_samples)),
  function(i) {
    rlang::try_fetch(
      {
        result = Screen(
          matched_bulk = bulk,
          sc_data = seurat,
          phenotype = pheno,
          label_type = "grade",
          phenotype_class = "binary",
          screen_method = "PIPET",
          distance = arg_samples$distance[i], # select_alpha will be used
          nPerm = as.integer(arg_samples$nPerm[i]),
          log2FC = arg_samples$log2FC[i]
        )

        data = data.frame(
          pos_cell = (result$scRNA_data$PIPET == "Positive")
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
  file = "pipet_random_search.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("pipet random search completed."))

# * EOC细胞标记为肿瘤细胞，作为判断依据，注意有肿瘤细胞转移的情况
pipet_random_search = data.table::fread(
  "pipet_random_search.csv",
)
pipet_random_search$benchmark = setNames(
  grepl("EOC", seurat$cell_subtype),
  colnames(seurat)
)
# * 计算指标
source(file.path(usethis::proj_path(), "Tmp/cross_args_ite/ComputeMetrics.R"))
metrics = ComputeMetrics(pipet_random_search)
t_metrics <- data.table::transpose(metrics)
colnames(t_metrics) <- rownames(metrics)
arg_samples = cbind(arg_samples, t_metrics)
data.table::fwrite(
  arg_samples,
  file = "pipet_arg_samples.csv",
  row.names = TRUE
)

# * 图
p <- ggplot(
  arg_samples,
  aes(
    x = nPerm,
    y = log2FC,
    fill = F1,
    shape = distance
  )
) +
  ggplot2::facet_wrap(~distance) +
  ggbeeswarm::geom_quasirandom(
    size = 6,
    alpha = 0.9,
    color = "black",
    method = "quasirandom", # 或 "swarm", "quasirandom"
    groupOnX = TRUE # 按x轴分组躲避
  ) +
  scale_shape_manual(values = c(21, 21, 21, 21, 21, 21)) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "F1"
  ) +
  labs(
    title = "Validation of the Screening Efficiency of degas under Random Parameters",
    subtitle = "x = nPerm, y = log2FC",
    x = "nPerm",
    y = "log2FC"
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
  filename = "pipet_acc.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)
