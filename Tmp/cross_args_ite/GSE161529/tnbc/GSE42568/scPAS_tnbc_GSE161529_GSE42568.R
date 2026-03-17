library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(BiocParallel)

# ! BRCA
# ! sc - GSE161529
# ! bulk - GSE42568
# ! survival

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/GSE161529/her2/GSE42568"
)


data_path = "/home/data/sigbridger/benchmark_data/brca"
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
  file.path(data_path, "seurat_tnbc.qs"),
  nthreads = 4
)

bulk = qs::qread(
  file.path(data_path, "brca_bulkdata_GSE42568.qs")
)
!anyNA(bulk)

pheno = qs::qread(file.path(data_path, "brca_pheno_GSE42568.qs"))

surv_data = select(
  pheno,
  "overall survival time_days:ch1",
  "overall survival event:ch1"
  # ,    "geo_accession"
) %>%
  rename("time" = 1, "status" = 2) %>%
  # tibble::column_to_rownames(var = "geo_accession") %>%
  filter(time != "NA" & status != "NA") %>%
  # mutate(
  #     status = case_when(
  #         status == "Alive" ~ 0,
  #         status == "Death" ~ 1
  #     )
  # ) %>%
  mutate_all(~ as.numeric(.))

# check
bulk = bulk[, rownames(surv_data)]
all(colnames(bulk) == rownames(surv_data))


# ! -------------------- scPAS ---------------------

# * random search, 100 times
set.seed(42568 + 1)
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
        scpas_result = Screen(
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
rownames(all_results) = colnames(seurat)

data.table::fwrite(
  all_results,
  file = "scpas_random_search.csv",
  row.names = TRUE
)

vroom::vroom_write(arg_samples, "scpas_arg_samples.csv", delim = ",")

cli::cli_alert_success(crayon::green("scPAS random search completed."))

# PID = 3761737
