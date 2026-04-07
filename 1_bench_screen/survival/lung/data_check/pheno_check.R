data_path <- "/home/data/sigbridger/benchmark_data/lung"

GSE3141_surv <- qs::qread(file.path(data_path, "GSE3141_surv_pheno.qs"))

GSE31210_surv <- qs::qread(
  file = file.path(data_path, "GSE31210_surv_pheno.qs")
)

GSE8894_surv <- qs::qread(
  file = file.path(data_path, "GSE8894_surv_pheno.qs")
)
TCGA_LUAD_surv <- qs::qread(
  file = file.path(data_path, "TCGA_LUAD_surv_pheno.qs")
)
