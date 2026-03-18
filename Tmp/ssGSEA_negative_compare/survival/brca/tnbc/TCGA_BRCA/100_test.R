setwd(usethis::proj_path())

scores <- qs::qread(
  "Tmp/ssGSEA_negative_compare/brca/tnbc_Sample_100_ssgsea_score.qs",
  nthreads = 4L
)


