library(dplyr)
library(data.table)
library(ScPP)
library(BiocParallel)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/survival/brca/her2/GSE42568"
)

data_path = "/home/data/sigbridger/benchmark_data/brca"
devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "seurat_her2.qs"))

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


gene_list <- marker_Survival2(bulk, surv_data)

# * random search, 51 times
set.seed(42568)
probs = round(seq(0.01, 0.5, length.out = 49), 2)

# ! estimate_cutoff and Log2FC_cutoff 不可在Survival表型中使用

# ! 并行设置
param <- MulticoreParam(workers = 4L)
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
# 2765512
