library(dplyr)
library(Seurat)
# library(BiocParallel)
library(scAB)
# library(SigBridgeR)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/binary/ov/GSE9891"
)

# devtools::document('/home/yyx/R/Project/R_code/SigBridgeR')
# data_path = "/home/data/sigbridger/benchmark_data/ov"

# # * load data
# seurat = qs::qread(file.path(data_path, "hgsoc_GSE165897_seurat.qs"))

# bulk = qs::qread(
#     file.path(data_path, "ov_bulkdata_GSE9891.qs"),
#     nthreads = 4
# )

# pheno = qs::qread(file.path(data_path, "ov_pheno_GSE9891.qs"))

# pheno_bi = setNames(
#     case_when(
#         pheno$characteristics_ch1.1 == "Type : LMP" ~ 0,
#         pheno$characteristics_ch1.1 == "Type : Malignant" ~ 1
#     ),
#     pheno$geo_accession
# )

# # ! 匹配
# bulk = bulk[, names(pheno_bi)]

# if (any(colnames(bulk) != names(pheno_bi))) {
#     stop("bulk and pheno_bi not match")
# }
# if (anyNA(bulk)) {
#     stop("bulk has NA")
# }
# if (anyNA(pheno_bi)) {
#     stop("pheno_bi has NA")
# }

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

# save_path <- "/home/data/sigbridger/method_compare/binary/ov/GSE9891"
# if (!dir.exists(save_path)) {
#     dir.create(save_path)
# }

# qs::qsave(scab_obj, file.path(save_path, "scab_obj.qs"))

scab_obj <- qs::qread(
  '/home/data/sigbridger/method_compare/binary/ov/GSE9891/scab_obj.qs'
)

para_list <- select_alpha.optimized(
  Object = scab_obj,
  method = 'binary',
  K = 2,
  cross_k = 5,
  para_1_list = alpha_samples,
  para_2_list = alpha_samples,
  parallel = F,
  workers = 2,
  verbose = TRUE,
  seed = 123
)

alpha <- para_list$para$alpha_1
alpha_2 <- para_list$para$alpha_2

cli::cli_alert_success(crayon::green(
  "Optimal parameters: alpha = {alpha}, alpha_2 = {alpha_2}"
))

qs::qsave(para_list, 'scAB_para_list.qs')
# 1853434
