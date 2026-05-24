# ! TCGA_LUAD
library(dplyr)

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
# * Screen

# * random search, 50 times
set.seed(123)
arg_samples <- data.frame(
  "nfeature" = sample(seq(500, 5000, by = 100), 50, replace = TRUE), # 第1维
  "imputation" = sample(c("None", "KNN", "ALRA"), 50, replace = TRUE), # 第2维
  "independent" = sample(c(FALSE, TRUE), 50, replace = TRUE) # 第3维
) %>%
  add_row(nfeature = 2000, imputation = "None", independent = TRUE)


# * run scpas
res_list <- lapply(
  seq_len(nrow(arg_samples)),
  function(i) {
    nfeature_i = arg_samples[i, "nfeature"][[1]]
    imputation_i = arg_samples[i, "imputation"][[1]]
    independent_i = arg_samples[i, "independent"][[1]]

    if (imputation_i == "None") {
      scpas_result = SigBridgeR::Screen(
        matched_bulk = bulk,
        sc_data = seurat,
        phenotype = pheno_bi,
        label_type = glue::glue("Tumor"),
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
      scpas_result = SigBridgeR::Screen(
        matched_bulk = bulk,
        sc_data = seurat,
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
  file = "stats/scpas_label_mat.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("(1) scPAS random search completed."))

# ! TCGA_LUAD
