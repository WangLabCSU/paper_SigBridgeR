setwd(usethis::proj_path())

scores <- qs::qread(
  "Tmp/ssGSEA_negative_compare/luad/luad_Sample_100_ssgsea_score.qs",
  nthreads = 4L
)



seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/TCGA-LUAD/tcga_luad_merged_seurat.qs",
  nthreads = 4L
)
meta <- seurat[[]]