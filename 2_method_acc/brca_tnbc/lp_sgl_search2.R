# ! GSE42568

setwd(file.path(usethis::proj_path(), "2_method_acc/brca_tnbc"))
library(dplyr)

# * Load Data
data_dir <- "/home/data/sigbridger/benchmark_data/brca"

sc_data <- qs::qread(file.path(data_dir, "seurat_tnbc.qs"), nthreads = 8L)

bulk <- qs::qread(
  file.path(data_dir, "brca_bulkdata_GSE42568.qs"),
  nthreads = 2L
)

cli::cli_alert_info("bulk data loaded: dim = ({.val {dim(bulk)}})")

pheno <- qs::qread(file.path(data_dir, "brca_pheno_GSE42568.qs"))

cm_samples <- intersect(rownames(pheno), colnames(bulk))


pheno_bi <- setNames(
  ifelse(pheno$`tissue:ch1` == "breast cancer", 1L, 0L),
  cm_samples
)
bulk <- bulk[, names(pheno_bi)]


cli::cli_alert_info("pheno data loaded: 1~tumor, 0~normal")
table(pheno_bi)

if (!all(names(pheno_bi) == colnames(bulk))) {
  stop("pheno_bi and bulk not match")
}

seurat_tumor <- readRDS(
  "/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_TNBCTum.rds"
)
tumor_cells <- rownames(seurat_tumor@meta.data)
benchmark_label <- colnames(sc_data) %in% tumor_cells

# * Screen

source("../rlogunif.R")

set.seed(123)
arg_samples <- data.frame(
  alpha = rlogunif(50, 0.005), # 第1维
  resolution = sample(seq(0, 1, 0.01), 50, replace = TRUE),
  nfold = sample(seq(5, 20, 1), 50, replace = TRUE)
) %>%
  dplyr::add_row(alpha = 0.5, resolution = 0.6, nfold = 5) # default parameters


# * run LP_SGL with error handling
res_list <- lapply(
  seq_len(nrow(arg_samples)),
  function(i) {
    cli::cli_h1("{i} / {nrow(arg_samples)}")
    result <- SigBridgeR::Screen(
      matched_bulk = bulk,
      sc_data = sc_data,
      phenotype = pheno_bi,
      label_type = glue::glue("process_{i}"),
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

    data <- data.frame(
      pos_cell = (result$scRNA_data$LP_SGL == "Positive")
    )
    colnames(data) = glue::glue("process_{i}")

    gc(verbose = FALSE)

    # 返回包含索引和结果的数据框
    return(data)
  }
)


# 合并所有结果
all_results <- do.call(cbind, res_list)
rownames(all_results) = colnames(seurat) # each cell is a row

data.table::fwrite(
  all_results,
  file = "stats/lp_sgl_label_mat1.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("(1)lpsgl random search completed."))

# ! GSE42568
