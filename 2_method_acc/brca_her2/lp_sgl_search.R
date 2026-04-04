setwd(file.path(usethis::proj_path(), "2_method_acc/brca_her2"))


# * Load Data
data_dir <- "/home/data/sigbridger/benchmark_data/brca"

sc_data <- qs::qread(file.path(data_dir, "seurat_her2.qs"), nthreads = 8L)

bulk <- qs::qread(file.path(data_dir, "brca_bulkdata_TCGA.qs"), nthreads = 2L)
bulk <- log2(bulk + 1)
cli::cli_alert_info("bulk data loaded: dim = ({.val {dim(bulk)}})")

pheno <- qs::qread(file.path(data_dir, "brca_pheno_TCGA.qs"))

cm_samples <- intersect(pheno$sample, colnames(bulk))

pheno_bi <- pheno %>%
  mutate(sample_type = substr(sample, 14, 15)) %>%
  filter(sample_type %in% c("01", "11"), sample %in% cm_samples) %>%
  mutate(sample_type = as.integer(sample_type == "01")) %>%
  {
    setNames(.$sample_type, .$sample)
  }

bulk <- bulk[, names(pheno_bi)]


cli::cli_alert_info("pheno data loaded: 1~tumor, 0~normal")
table(pheno_bi)

if (!all(names(pheno_bi) == colnames(bulk))) {
  stop("pheno_bi and bulk not match")
}

seurat_tumor <- readRDS(
  "/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_HER2Tum.rds"
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


# * run scab with error handling
res_list <- lapply(
  seq_len(nrow(arg_samples)),
  function(i) {
    result = Screen(
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

    data = data.frame(
      pos_cell = (result$scRNA_data$LP_SGL == "Positive")
    )
    colnames(data) = glue::glue("process_{i}")

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
