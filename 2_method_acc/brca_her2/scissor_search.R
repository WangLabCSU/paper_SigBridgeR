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

# ? warmup

tmp <- SigBridgeR::Screen(
  matched_bulk = bulk,
  sc_data = sc_data,
  phenotype = pheno_bi,
  label_type = "tumor",
  phenotype_class = "binary",
  screen_method = "Scissor",
  alpha = 0.9,
  cutoff = 0.2,
  path2save_scissor_inputs = "TCGA_BRCA_her2_scissor_cache.RData"
)
rm(tmp)


alpha <- c(0.001, seq(0.05, 0.95, 0.05))
cutoff <- seq(0.05, 0.5, 0.05)

results <- lapply(cutoff, \(c) {
  res <- SigBridgeR::Screen(
    matched_bulk = bulk,
    sc_data = sc_data,
    phenotype = pheno_bi,
    label_type = glue::glue("process_{c}"),
    phenotype_class = "binary",
    screen_method = "scissor",
    alpha = alpha,
    cutoff = c,
    path2load_scissor_cache = "TCGA_BRCA_her2_scissor_cache.RData"
  )

  pos_ratio = (res$scRNA_data$scissor == "Positive")
  pos = data.frame(pos = pos_ratio)
  colnames(pos) <- glue::glue("process_{c}")
  pos
})

results <- dplyr::bind_cols(results)
rownames(results) = colnames(sc_data)


data.table::fwrite(
  res,
  file = "stats/scissor_label_mat.csv",
  row.names = TRUE
)

cli::cli_alert_info("scissor_label_mat saved")

# ! TCGA_BRCA
