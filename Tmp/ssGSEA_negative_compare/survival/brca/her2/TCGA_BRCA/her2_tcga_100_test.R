setwd(usethis::proj_path())

scores <- qs::qread(
  "Tmp/ssGSEA_negative_compare/brca/her2_Sample_100_ssgsea_score.qs",
  nthreads = 4L
)

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_data/brca/HER2/tcga_her2_merged_seurat.qs",
  nthreads = 4L
)

meta <- seurat[[]]
