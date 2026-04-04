data_dir <- "/home/data/sigbridger/benchmark_data/ov"

GSE32062_bulk <- qs::qread(file.path(
  data_dir,
  "ov_bulkdata_GSE32062_GPL6480.qs"
))
GSE32062_pheno <- qs::qread(file.path(
  data_dir,
  "ov_pheno_GSE32062.qs"
))
GSE32062_surv <- GSE32062_pheno %>%
  tibble::column_to_rownames("Sample_ID") %>%
  dplyr::select(`OS (M)`, `Death (1)`) %>%
  dplyr::rename(time = 1, status = 2)

GSE140082_bulk <- qs::qread(file.path(
  data_dir,
  "ov_bulkdata_GSE140082.qs"
))
GSE140082_pheno <- qs::qread(file.path(data_dir, "ov_pheno_GSE140082.qs"))
GSE140082_surv <- GSE140082_pheno %>%
  dplyr::select(`final_ostm:ch1`, `final_osid:ch1`) %>%
  dplyr::rename(time = 1, status = 2) %>%
  dplyr::mutate_all(as.integer)
