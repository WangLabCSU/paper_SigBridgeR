setwd(file.path(usethis::proj_path(), "2_method_acc/brca_her2"))


# * Load Data
data_dir <- "/home/data/sigbridger/benchmark_data/brca"

sc_data <- qs::qread(file.path(data_dir, "seurat_her2.qs"), nthreads = 8L)

bulk <- qs::qread(file.path(data_dir, "brca_bulkdata_TCGA.qs"), nthreads = 2L)

cli::cli_alert_info("bulk data loaded: dim = ({.val {dim(bulk)}})")

pheno <- qs::qread(file.path(data_dir, "brca_pheno_TCGA.qs"))

pheno_samples <- rownames(pheno)

cm_samples <- intersect(pheno_samples, colnames(bulk))
bulk <- bulk[, cm_samples]
pheno_bi <- ifelse(pheno[cm_samples, ]$sample %>% endsWith("01"), 1, 0)
names(pheno_bi) <- cm_samples

cli::cli_alert_info("pheno data loaded: 1~tumor, 0~normal")
table(pheno_bi)

if (!all(names(pheno_bi) == colnames(bulk))) {
  stop("pheno_bi and bulk not match")
}

# * Screen
