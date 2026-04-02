library(dplyr)
library(data.table)
library(ScPP)
# library(furrr)
library(BiocParallel)


# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/binary/brca/her2/GSE162228"
)

data_path <- "/home/data/sigbridger/benchmark_data/brca"
devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat <- qs::qread(file.path(data_path, "seurat_her2.qs"))

bulk <- qs::qread(
  file.path(data_path, "brca_bulkdata_GSE162228.qs")
)

pheno <- qs::qread(file.path(
  data_path,
  "brca_pheno_GSE162228.qs"
))

pheno_bi <- setNames(
  case_when(
    pheno$`relapse status:ch1` == "relapse" ~ 1,
    pheno$`relapse status:ch1` == "non-relapse" ~ 0
  ),
  pheno$geo_accession
)


if (!all(colnames(bulk) == names(pheno_bi))) {
  stop("bulk and pheno_bi not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}
if (anyNA(pheno_bi)) {
  stop("pheno_bi has NA")
}

# --------------------------------------------------------------------------------------------
# ! scPP random search, 50 times

# * random search, 51 times
set.seed(12345)
probs_sample <- round(seq(0.01, 0.5, length.out = 49), 2)

arg_sample <- data.frame(
  prob = sample(probs_sample, 50, replace = TRUE),
  Log2FC_cutoff = round(runif(50), 3)
) %>%
  add_row(prob = 0.2, Log2FC_cutoff = 0.585) # default


# ------------------------------------------------------------------------------------------

# ! 并行设置
# future::plan("multisession", workers = 3L)
# options(future.globals.maxSize = 1024^3 * 5L)

# param <- BiocParallel::MulticoreParam(workers = 2L)
# BiocParallel::register(param) # 注册为默认后端

res_list <- purrr::map(
  seq_len(nrow(arg_sample)),
  function(i) {
    prob_i <- arg_sample$prob[i]
    Log2FC_cutoff_i <- arg_sample$Log2FC_cutoff[i]

    scpp_result = Screen(
      bulk,
      seurat,
      pheno_bi,
      screen_method = "scPP",
      label_type = "relapse",
      phenotype_class = "binary",
      ref_group = 0L,
      Log2FC_cutoff = Log2FC_cutoff_i,
      probs = prob_i,
      parallel = FALSE
    )

    pos_cell <- (scpp_result$scPP == "Positive")

    data = data.frame(
      pos_cell = pos_cell
    )
    if (nrow(data) < ncol(seurat)) {
      cli::cli_alert_warning(c(
        " " = "scPP result is not complete, maybe not suitable for this parameter pair, using default FALSE.",
        ">" = "prob = {prob_i}, Log2FC_cutoff = {Log2FC_cutoff_i}, process_{i}"
      ))
      data = data.frame(pos_cell = rep(FALSE, ncol(seurat)))
    }
    colnames(data) = glue::glue(
      "process_{prob_i}_{Log2FC_cutoff_i}"
    )
    gc(verbose = FALSE)

    # 返回包含索引和结果的数据框
    return(data)
  }
  # ,
  # .options = furrr::furrr_options(
  #     packages = c("SigBridgeR", "glue", "ScPP"),
  #     seed = TRUE,
  #     globals = TRUE
  # )
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

# 617004
