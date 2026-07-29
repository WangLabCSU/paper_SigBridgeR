# ! run this after `0_bench_data_preprocess/brca_bulk_data.R`
# GSE32062, GSE140082, GSE9891
setwd(usethis::proj_path())


data_dir <- "/home/data/sigbridger/benchmark_data/ov"
gse8894 <- qs::qread(file.path(data_dir, "ov_bulkdata_GSE9891.qs"))
gse32062 <- qs::qread(file.path(data_dir, "ov_bulkdata_GSE32062_GPL6480.qs")) # normalized
gse140082 <- qs::qread(file.path(data_dir, "ov_bulkdata_GSE140082.qs")) # TPM
