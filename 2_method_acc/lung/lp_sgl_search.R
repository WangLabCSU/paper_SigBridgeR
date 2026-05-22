# ! TCGA_LUAD

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  file.path(usethis::proj_path(), "2_method_acc/lung")
)
data_path <- "/home/data/sigbridger/benchmark_data/lung"

# * load data
seurat <- qs::qread(file.path(data_path, "luad_GSE123902_seurat.qs"))

bulk <- qs::qread(
  file.path(data_path, "TCGA_LUAD_bulkdata.qs")
)
bulk <- log2(bulk + 1)

pheno <- qs::qread(file.path(data_path, "TCGA_LUAD_pheno.qs"))

pheno_bi <- mutate(pheno, sample_type = substr(pheno$sample, 14, 15)) %>%
  select(sample, sample_type) %>%
  filter(sample_type %in% c("01", "11")) %>%
  mutate(sample_type = as.integer(sample_type == "01"))
pheno_bi <- setNames(pheno_bi$sample_type, pheno_bi$sample)

bulk <- bulk[, names(pheno_bi)]

if (!all(colnames(bulk) == names(pheno_bi))) {
  stop("bulk and pheno_bi not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}
if (anyNA(pheno_bi)) {
  stop("pheno_bi has NA")
}

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

# ! TCGA_LUAD
