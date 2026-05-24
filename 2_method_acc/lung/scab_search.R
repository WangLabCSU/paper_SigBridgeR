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

# future::plan(future::multicore, workers = 5L)
SigBridgeR::setThreads(4L)


# ! To avoid recomputing, file cache is used
if (!dir.exists("stats/scab1")) {
  dir.create("stats/scab1", recursive = TRUE)
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

if (!file.exists("stats/scab_label_mat1.csv")) {
  res_list <- lapply(
    seq_len(nrow(arg_samples)),
    function(i) {
      cli::cli_h1("{i} / {nrow(arg_samples)}")

      cache_save_path <- file.path(
        "stats/scab1",
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
        parallel = FALSE,
        verbose = TRUE
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
    all_results,
    file = "stats/scab_label_mat1.csv",
    row.names = TRUE
  )
}


cli::cli_alert_success(crayon::green("(1) scab random search completed."))

# ! TCGA_LUAD
