# ! run this after `0_bench_data_preprocess/brca_bulk_data.R`
# GSE3141, GSE8894, GSE31210, TCGA_LUAD
setwd(usethis::proj_path())


library(zeallot)
data_dir <- "/home/data/sigbridger/benchmark_data/lung"

c(gse3141, gse8894, gse31210, tcga_luad) %<-%
  purrr::map(
    c("GSE3141", "GSE8894", "GSE31210", "TCGA_LUAD"),
    \(x) {
      if (!x == "TCGA_LUAD") {
        qs::qread(file.path(data_dir, glue::glue("lung_bulkdata_{x}.qs")))
      } else {
        qs::qread(file.path(data_dir, "TCGA_LUAD_bulkdata.qs"))
      }
    }
  )

# gse3141 - TPM
# gse8894 - log2 GC-RMA signal intensity
# gse31210 - TPM

source("0_bench_data_preprocess/Counts2TPM.R")
gene_length_cmb <- data.table::fread(
  "0_bench_data_preprocess/gencode_v33_tcga_gene_length.csv"
)
gene_length_list <- setNames(
  gene_length_cmb$gene_length,
  gene_length_cmb$gene_id
)
tcga_luad_tpm <- Counts2TPM(as.matrix(tcga_luad), gene_length_list)
qs::qsave(
  tcga_luad_tpm,
  file.path(data_dir, "TCGA_LUAD_bulkdata_tpm.qs"),
  nthreads = 4L
)
