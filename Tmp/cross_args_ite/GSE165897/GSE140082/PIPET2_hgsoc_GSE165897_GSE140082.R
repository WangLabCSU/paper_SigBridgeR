library(dplyr)
library(Seurat)
library(BiocParallel)
library(SigBridgeR)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/GSE165897/GSE140082"
)

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


# * run pipet with error handling
res_list <- lapply(
  seq_len(nrow(arg_samples)),
  function(i) {
    rlang::try_fetch(
      {
        result = Screen(
          matched_bulk = bulk,
          sc_data = seurat,
          phenotype = pheno,
          label_type = glue::glue("OS (M)_survival_{i}"),
          phenotype_class = "binary",
          screen_method = "PIPET",
          distance = arg_samples$distance[i], # select_alpha will be used
          nPerm = as.integer(arg_samples$nPerm[i]),
          log2FC = arg_samples$log2FC[i]
        )

        # qs::qsave(
        #     pipet_result,
        #     file = glue::glue("pipet_results/pipet_result_{i}.qs")
        # )

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

# -----------------------------------------------------------------------------------------------
# ! viz

pipet_random_search = data.table::fread(
  "pipet_random_search.csv",
)
pipet_random_search$benchmark = setNames(
  grepl("EOC", seurat$cell_subtype),
  colnames(seurat)
)

# * 计算指标
source(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/ComputeMetrics.R"
)
metrics = ComputeMetrics(pipet_random_search)
arg_samples$Accuracy = data.table::transpose(metrics)[[6]]

data.table::fwrite(
  arg_samples,
  file = "pipet_arg_samples.csv",
  row.names = TRUE
)
library(ggplot2)

p <- ggplot(
  arg_samples,
  aes(
    x = nPerm,
    y = log2FC,
    fill = Accuracy
    # 移除 shape 映射
  )
) +
  geom_point(
    size = 5,
    alpha = 0.9,
    color = "black",
    shape = 21, # 统一使用支持 fill 的形状
    position = position_jitter(width = 0.2, height = 0.15, seed = 42)
  ) +
  scale_fill_gradient(low = "white", high = "red", name = "Accuracy") +
  facet_wrap(~distance, ncol = 3) + # 按 arch 分面，6 个类别自动排列
  labs(
    title = "Validation of the Screening Efficiency of degas under Random Parameters",
    subtitle = "x = bag_depth, y = ff_depth",
    x = "bag_depth",
    y = "ff_depth"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    strip.text = element_text(size = 13, face = "bold") # 分面标签样式
  )

ggsave(
  filename = "pipet_random_search.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)
