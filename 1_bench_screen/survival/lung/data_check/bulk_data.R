data_dir <- "/home/data/sigbridger/benchmark_data/lung"

TCGA_LUAD_bulk <- qs::qread(file.path(data_dir, "TCGA_LUAD_bulkdata.qs")) # counts

GSE3141_bulk <- qs::qread(file.path(data_dir, "lung_bulkdata_GSE3141.qs")) # counts

GSE8894_bulk <- qs::qread(file.path(data_dir, "lung_bulkdata_GSE8894.qs"))

GSE31210_bulk <- qs::qread(file.path(data_dir, "lung_bulkdata_GSE31210.qs")) # counts

GSE3141_surv <- qs::qread(file.path(data_dir, "GSE3141_surv_pheno.qs"))
