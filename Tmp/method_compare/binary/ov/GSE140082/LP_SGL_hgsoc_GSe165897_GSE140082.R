library(dplyr)
library(Seurat)
library(BiocParallel)
library(SigBridgeR)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/GSE165897/GSE140082"
)

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


# * random search, 100 times
set.seed(123)
arg_samples <- data.frame(
  alpha = rlogunif(50, 0.005), # 第1维
  resolution = sample(seq(0, 1, 0.01), 50, replace = TRUE),
  nfold = sample(seq(5, 20, 1), 50, replace = TRUE)
) %>%
  add_row(alpha = 0.5, resolution = 0.6, nfold = 5) # default parameters


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
          label_type = glue::glue("grade_{i}"),
          phenotype_class = "binary",
          screen_method = "LP_SGL",
          alpha = arg_samples$alpha[i], # select_alpha will be used
          resolution = arg_samples$resolution[i],
          nfold = as.integer(arg_samples$nfold[i])
        )

        # qs::qsave(
        #     scab_result,
        #     file = glue::glue("scAB_results/scab_result_{i}.qs")
        # )

        data = data.frame(
          pos_cell = (result$scRNA_data$LP_SGL == "Positive")
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
  file = "lpsgl_random_search.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("lpsgl random search completed."))
