# ! TCGA_BRCA

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

set.seed(123)
alpha_samples <- Reduce(`*`, rep(c(2, 5), 2), init = 5e-4, accumulate = TRUE)
tred <- 1:10

save_path = '/home/data/sigbridger/method_compare/binary/brca/her2/TCGA_BRCA'

scAB_obj <- scAB::create_scAB.v5(
  Object = sc_data,
  bulk_dataset = bulk,
  phenotype = pheno_bi,
  method = "binary"
)

k <- scAB::select_K.optimized(
  Object = scAB_obj,
  K_max = 20L,
  repeat_times = 10L,
  maxiter = 2000L, # default in scAB
  seed = seed,
  verbose = verbose
)

scab_res = scAB::scAB.optimized(
  Object = scab_obj,
  K = k,
  alpha = alpha_samples,
  alpha_2 = alpha_samples
)

combind_res <- data.frame()
for (i in tred) {
  seurat_screened <- scAB::findSubset.optimized(
    Object = sc_data,
    scAB_Object = scab_res,
    tred = i
  )

  label <- data.frame(
    seurat_screened$scAB == "Positive"
  )
  colnames(label) <- paste("process", i, sep = "_")
  combind_res <- cbind(combind_res, label)
  gc()
}

rownames(combind_res) <- colnames(sc_data)

data.table::fwrite(
  all_results,
  file = "stats/scab_label_mat1.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("(1)scab random search completed."))
