# ! GSE42568

setwd(file.path(usethis::proj_path(), "2_method_acc/brca_her2"))
library(dplyr)


# * Load Data
data_dir <- "/home/data/sigbridger/benchmark_data/brca"

sc_data <- qs::qread(file.path(data_dir, "seurat_her2.qs"), nthreads = 8L)

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
  "/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_HER2Tum.rds"
)
tumor_cells <- rownames(seurat_tumor@meta.data)
benchmark_label <- colnames(sc_data) %in% tumor_cells

# * Screen

future::plan(future::multicore, workers = 5L)
# ! To avoid recomputing, file cache is used
if (!dir.exists("stats/scab2")) {
  dir.create("stats/scab2", recursive = TRUE)
}


set.seed(123)
alpha_samples <- Reduce(`*`, rep(c(2, 5), 2), init = 5e-4, accumulate = TRUE)

arg_samples <- data.frame(
  tred = sample(seq(0.5, 5, 0.5), 50, replace = TRUE),
  repeat_times = sample(5:20, 50, replace = TRUE)
) %>%
  dplyr::add_row(tred = 2L, repeat_times = 10L)

scAB_obj <- scAB::create_scAB.v5(
  Object = sc_data,
  bulk_dataset = bulk,
  phenotype = pheno_bi,
  method = "binary"
)

if (!file.exists("stats/scab_label_mat2.csv")) {
  res_list <- lapply(
    seq_len(nrow(arg_samples)),
    function(i) {
      cli::cli_h1("{i} / {nrow(arg_samples)}")

      cache_save_path <- file.path(
        "stats/scab2",
        glue::glue("process_{i}.csv")
      )
      if (file.exists(cache_save_path)) {
        cli::cli_alert("cache found, loading...")
        cache <- data.table::fread(cache_save_path)
        return(cache)
      }

      repeat_times_i <- arg_samples[i, "repeat_times"][[1]]
      tred_i <- arg_samples[i, "tred"][[1]]

      k <- scAB::select_K.optimized(
        Object = scAB_obj,
        K_max = 20L,
        repeat_times = repeat_times_i,
        maxiter = 2000L # default in scAB
      )

      para_list <- scAB::select_alpha.optimized(
        Object = scAB_obj,
        method = "binary",
        K = k,
        cross_k = 5,
        para_1_list = alpha_samples %||% c(0.01, 0.005, 0.001),
        para_2_list = alpha_samples %||% c(0.01, 0.005, 0.001),
        parallel = TRUE,
        verbose = FALSE
      )

      alpha <- para_list$para$alpha_1
      alpha_2 <- para_list$para$alpha_2

      scab_res <- scAB::scAB.optimized(
        Object = scAB_obj,
        K = k,
        alpha = alpha,
        alpha_2 = alpha_2
      )

      seurat_screened <- scAB::findSubset.optimized(
        Object = sc_data,
        scAB_Object = scab_res,
        tred = tred_i
      )

      label <- data.frame(
        seurat_screened$scAB == "Positive"
      )
      colnames(label) <- paste("process", i, sep = "_")

      # ! save cache
      data.table::fwrite(label, cache_save_path)

      gc()
      label
    }
  )

  gc()
  all_results <- do.call(cbind, res_list)
  rownames(all_results) = colnames(sc_data)

  data.table::fwrite(
    res_list,
    file = "stats/scab_label_mat2.csv",
    row.names = TRUE
  )
}

cli::cli_alert_success(crayon::green("(2) scab random search completed."))
