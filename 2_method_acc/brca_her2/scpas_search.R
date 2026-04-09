setwd(file.path(usethis::proj_path(), "2_method_acc/brca_her2"))
library(dplyr)


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
      scpas_result = Screen(
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
      scpas_result = Screen(
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
rownames(all_results) = colnames(sc_data)

data.table::fwrite(
  all_results,
  file = "stats/scpas_label_mat.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("(1) scPAS random search completed."))

# ! TCGA_BRCA
