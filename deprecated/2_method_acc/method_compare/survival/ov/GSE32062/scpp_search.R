library(dplyr)
library(data.table)
library(ScPP)
library(BiocParallel)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/survival/ov/GSE140082"
)

data_path = "/home/data/sigbridger/benchmark_data/ov"
devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "hgsoc_GSE165897_seurat.qs"))

bulk = qs::qread(
  file.path(data_path, "ov_bulkdata_GSE32062_GPL6480.qs")
)

pheno = qs::qread(file.path(data_path, "ov_pheno_GSE32062.qs"))

surv_data = select(
  pheno,
  "OS (M)",
  "Death (1)",
  "Sample_ID"
) %>%
  rename("time" = 1, "status" = 2) %>%
  tibble::column_to_rownames(var = "Sample_ID") %>%
  mutate_all(~ as.numeric(.))

if (!all(colnames(bulk) == rownames(surv_data))) {
  stop("bulk and surv_data not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}

# ? _____________________ scPP _____________________

gene_list <- marker_Survival2(bulk, surv_data)

# * random search, 51 times
set.seed(123)
probs = round(seq(0.01, 0.5, length.out = 49), 2)

# ! estimate_cutoff and Log2FC_cutoff 不可在Survival表型中使用

# ! 并行设置
param <- MulticoreParam(workers = 2)
register(param) # 注册为默认后端

res_list <- bplapply(
  probs,
  function(prob_i) {
    tryCatch(
      {
        scpp_result = FixedProbMode(
          sc_dataset = seurat,
          geneList = gene_list,
          prob = prob_i
        )

        pos_cell = (scpp_result$metadata$scPP == "Positive")

        data = data.frame(
          pos_cell = pos_cell
        )
        colnames(data) = glue::glue("process_{prob_i}")
        gc()

        # 返回包含索引和结果的数据框
        return(data)
      },
      error = function(e) {
        cli::cli_alert_danger(e$message)
        data = data.frame(
          pos_cell = FALSE
        )
        colnames(data) = glue::glue("process_{prob_i}")
        return(data)
      }
    )
  }
)


# *visualize
gc()
all_results <- do.call(cbind, res_list)
rownames(all_results) = colnames(seurat)

data.table::fwrite(
  all_results,
  file = "scpp_random_search.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("scPP random search completed."))
# 2964209
