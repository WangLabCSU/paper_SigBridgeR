# ! TCGA_LUAD

library(dplyr)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  file.path(usethis::proj_path(), "2_method_acc/lung")
)
data_path <- "/home/data/sigbridger/benchmark_data/lung"

# * load data
sc_data <- qs::qread(file.path(data_path, "luad_GSE123902_seurat.qs"))

bulk <- qs::qread(
  file.path(data_path, "TCGA_LUAD_bulkdata_tpm.qs")
)

pheno <- qs::qread(file.path(data_path, "TCGA_LUAD_pheno.qs"))

pheno_bi <- mutate(pheno, sample_type = substr(pheno$sample, 14, 15)) %>%
  select(sample, sample_type) %>%
  filter(sample_type %in% c("01", "11")) %>%
  mutate(sample_type = as.integer(sample_type == "01"))
pheno_bi <- setNames(pheno_bi$sample_type, pheno_bi$sample)

cm_samples <- intersect(colnames(bulk), names(pheno_bi))

bulk <- bulk[, cm_samples]
pheno_bi <- pheno_bi[cm_samples]

if (!all(colnames(bulk) == names(pheno_bi))) {
  stop("bulk and pheno_bi not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}
if (anyNA(pheno_bi)) {
  stop("pheno_bi has NA")
}

# ? warmup

if (!file.exists("TCGA_LUAD_tpm_lung_scissor_cache.RData")) {
  tmp <- SigBridgeR::Screen(
    matched_bulk = bulk,
    sc_data = sc_data,
    phenotype = pheno_bi,
    label_type = "tumor",
    phenotype_class = "binary",
    screen_method = "Scissor",
    alpha = 0.9,
    cutoff = 0.2,
    path2save_scissor_inputs = "TCGA_LUAD_tpm_lung_scissor_cache.RData"
  )
  rm(tmp)
}

alpha <- c(0.001, seq(0.05, 0.95, 0.05))
cutoff <- seq(0.05, 0.5, 0.05)

results <- lapply(cutoff, \(c) {
  res <- SigBridgeR::Screen(
    matched_bulk = bulk,
    sc_data = sc_data,
    phenotype = pheno_bi,
    label_type = glue::glue("process_{c}"),
    phenotype_class = "binary",
    screen_method = "Scissor",
    alpha = alpha,
    cutoff = c,
    path2load_scissor_cache = "TCGA_LUAD_tpm_lung_scissor_cache.RData"
  )

  pos_ratio = (res$scRNA_data$scissor == "Positive")
  pos = data.frame(pos = pos_ratio)
  colnames(pos) <- glue::glue("process_{c}")
  pos
})

results <- dplyr::bind_cols(results)
rownames(results) = colnames(sc_data)


data.table::fwrite(
  results,
  file = "stats/scissor_label_mat1.csv",
  row.names = TRUE
)

cli::cli_alert_info("scissor_label_mat saved")

# ! TCGA_LUAD
