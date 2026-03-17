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

# * random search, 100 times
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
for (j in c(1, 2)) {
  # j = 2
  sub_arg_sample <- split_arg_sample[[j]]

  res_list <- lapply(
    seq_len(nrow(sub_arg_sample)),
    function(i) {
      if (j == 1 && i == 10 || j == 2 && i == 3) {
        return(data.frame(pos_cell = rep(FALSE, ncol(seurat))))
      }

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
      colnames(data) = glue::glue("process_{j}{i}")

      data.table::fwrite(
        data,
        file = file.path(
          save_path,
          paste0(j, i, "_scpas_random_search.csv")
        ),
        row.names = TRUE
      )
      gc()
      cli::cli_alert_success("Finished j={j}, i={i} job")

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

#
