# ! GSE42568

setwd(file.path(usethis::proj_path(), "2_method_acc/brca_tnbc"))


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


bulk <- bulk[, names(pheno_bi)]
pheno_bi <- setNames(
  ifelse(pheno$`tissue:ch1` == "breast cancer", 1L, 0L),
  cm_samples
)

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

set.seed(42568)
arg_samples <- data.frame(
  "nfeature" = sample(seq(500, 5000, by = 100), 50, replace = TRUE), # 第1维
  "imputation" = sample(c("None", "KNN", "ALRA"), 50, replace = TRUE), # 第2维
  "independent" = sample(c(FALSE, TRUE), 50, replace = TRUE) # 第3维
) %>%
  dplyr::add_row(nfeature = 2000, imputation = "None", independent = TRUE)


# * run scpas
res_list <- lapply(
  seq_len(nrow(arg_samples)),
  function(i) {
    nfeature_i = arg_samples[i, "nfeature"][[1]]
    imputation_i = arg_samples[i, "imputation"][[1]]
    independent_i = arg_samples[i, "independent"][[1]]

    if (imputation_i == "None") {
      scpas_result = Screen(
        matched_bulk = bulk,
        sc_data = sc_data,
        phenotype = pheno_bi,
        label_type = glue::glue("process_{i}"),
        phenotype_class = "binary",
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
        sc_data = sc_data,
        phenotype = pheno_bi,
        label_type = glue::glue("Tumor"),
        phenotype_class = "binary",
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
    colnames(data) = glue::glue("process_{i}")
    gc()

    # 返回包含索引和结果的数据框
    return(data)
  }
)


# *visualize
gc()
all_results <- do.call(cbind, res_list)
rownames(all_results) = colnames(seurat)

data.table::fwrite(
  all_results,
  file = "stats/scpas_label_mat2.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("(2)scpas random search completed."))

# ! GSE42568
