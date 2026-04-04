data_dir <- "/home/data/sigbridger/benchmark_data/lung"

GSE3141_pheno <- data.table::fread(file.path(
  data_dir,
  "lung.cancer.adeno.gse3141.hgu133plus2_entrezcdf.tsv"
))
GSE3141_surv <- GSE3141_pheno %>%
  tibble::column_to_rownames("Array") %>%
  dplyr::select(OS_Time, OS_Status) %>%
  dplyr::rename(time = 1, status = 2)

qs::qsave(GSE3141_surv, file.path(data_dir, "GSE3141_surv_pheno.qs"))

GSE8894_pheno <- data.table::fread(file.path(
  data_dir,
  "lung.cancer.adeno.gse8894.gpl570.tsv"
))
GSE8894_surv <- GSE8894_pheno %>%
  tibble::column_to_rownames("Array") %>%
  dplyr::select(OS_Time, OS_Status) %>%
  dplyr::rename(time = 1, status = 2)

qs::qsave(GSE8894_surv, file.path(data_dir, "GSE8894_surv_pheno.qs"))

GSE31210_pheno <- data.table::fread(file.path(
  data_dir,
  "lung.cancer.adeno.gse31210.hgu133plus2_entrezcdf.tsv"
))
GSE31210_surv <- GSE31210_pheno %>%
  tibble::column_to_rownames("Array") %>%
  dplyr::select(OS_Time, OS_Status) %>%
  dplyr::rename(time = 1, status = 2)

qs::qsave(GSE31210_surv, file.path(data_dir, "GSE31210_surv_pheno.qs"))
