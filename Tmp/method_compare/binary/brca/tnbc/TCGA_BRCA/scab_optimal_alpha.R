library(dplyr)
library(Seurat)
# library(BiocParallel)
library(scAB)
# library(SigBridgeR)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/binary/brca/tnbc/TCGA_BRCA"
)

# data_path = "/home/data/sigbridger/benchmark_data/brca"
# devtools::document("~/R/Project/R_code/SigBridgeR")

# # * load data
# seurat = qs::qread(file.path(data_path, "seurat_tnbc.qs"))

# bulk = qs::qread(
#     file.path(data_path, "brca_bulkdata_TCGA.qs")
# )

# pheno = qs::qread(
#     file.path(data_path, "brca_pheno_TCGA.qs"),
#     nthreads = 4
# )

# pheno_bi = mutate(
#     pheno,
#     sample_type = substr(pheno$sample, 14, 15)
# ) %>%
#     select(sample, sample_type) %>%
#     filter(sample_type %in% c("01", "11")) %>%
#     mutate(sample_type = ifelse(sample_type == "01", 1, 0))

# pheno_bi = setNames(pheno_bi$sample_type, pheno_bi$sample)

# cm_samples = intersect(names(pheno_bi), colnames(bulk))
# bulk = bulk[, cm_samples]
# pheno_bi = pheno_bi[cm_samples]

# if (!all(colnames(bulk) == names(pheno_bi))) {
#     stop("bulk and pheno_bi not match")
# }
# if (anyNA(bulk)) {
#     stop("bulk has NA")
# }
# if (anyNA(pheno_bi)) {
#     stop("pheno_bi has NA")
# }
# # ----------------------------scAB--------------------------------------------------------------------------------

set.seed(123)
alpha_samples <- Reduce(`*`, rep(c(2, 5), 2), init = 5e-4, accumulate = TRUE)

# scab_obj <- create_scAB.v5(
#     Object = seurat,
#     bulk_dataset = bulk,
#     phenotype = pheno_bi,
#     method = 'binary',
#     verbose = TRUE
# )

# k <- scAB::select_K.optimized(
#     Object = scab_obj,
#     K_max = 20L,
#     repeat_times = 10L,
#     maxiter = 2000L, # default in scAB
#     seed = 123,
#     verbose = TRUE
# )

# cli::cli_alert_success(crayon::green(
#     "Optimal K = {k}"
# ))

save_path <- "/home/data/sigbridger/method_compare/binary/brca/tnbc/TCGA_BRCA"
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# qs::qsave(scab_obj, file.path(save_path, "scab_obj.qs"))

scab_obj <- qs::qread(
  '/home/data/sigbridger/method_compare/binary/brca/tnbc/TCGA_BRCA/scab_obj.qs'
)

para_list <- select_alpha.optimized(
  Object = scab_obj,
  method = 'binary',
  K = 6L,
  cross_k = 5,
  para_1_list = alpha_samples,
  para_2_list = alpha_samples,
  parallel = FALSE,
  verbose = TRUE,
  seed = 123
)

alpha <- para_list$para$alpha_1
alpha_2 <- para_list$para$alpha_2

cli::cli_alert_success(crayon::green(
  "Optimal parameters: alpha = {alpha}, alpha_2 = {alpha_2}"
))

qs::qsave(para_list, 'scAB_para_list.qs')
# 878521
