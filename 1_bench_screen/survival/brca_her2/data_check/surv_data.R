data_path <- "/home/data/sigbridger/benchmark_data/brca"

GSE42568_pheno <- qs::qread(file.path(data_path, "brca_pheno_GSE42568.qs"))
GSE42568_bulk <- qs::qread(file.path(data_path, "brca_bulkdata_GSE42568.qs"))

surv_data <- GSE42568_pheno %>%
  dplyr::select(
    `overall survival time_days:ch1`,
    `overall survival event:ch1`
  ) %>%
  dplyr::filter(
    !is.na(`overall survival time_days:ch1`) &
      !is.na(`overall survival event:ch1`) &
      `overall survival time_days:ch1` != "NA"
  ) %>%
  dplyr::rename(time = 1, status = 2)


GSE162228_pheno <- qs::qread(file.path(data_path, "brca_pheno_GSE162228.qs"))
GSE162228_bulk <- qs::qread(file.path(data_path, "brca_bulkdata_GSE162228.qs"))

surv_data <- GSE162228_pheno %>%
  dplyr::select(`overall survival (years):ch1`, `alive:ch1`) %>%
  dplyr::rename(time = 1, status = 2) %>%
  dplyr::mutate(
    status = dplyr::case_when(status == "Alive" ~ 1L, status == "Death" ~ 0L)
  )
