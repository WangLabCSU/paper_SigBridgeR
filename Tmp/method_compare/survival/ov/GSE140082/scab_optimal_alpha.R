library(dplyr)
library(Seurat)
# library(BiocParallel)
library(scAB)
# library(SigBridgeR)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/survival/ov/GSE140082"
)

devtools::document('/home/yyx/R/Project/R_code/SigBridgeR')
data_path = "/home/data/sigbridger/benchmark_data/ov"

# * load data
seurat = qs::qread(file.path(data_path, "hgsoc_GSE165897_seurat.qs"))

bulk = qs::qread(
  file.path(data_path, "ov_bulkdata_GSE140082.qs")
)
# GSM4153781 has NA values, find median and replace NA
median = median(bulk[, "GSM4153781"], na.rm = TRUE)
na_indices <- which(is.na(bulk[, "GSM4153781"]))
bulk[na_indices, "GSM4153781"] <- median

pheno = qs::qread(file.path(data_path, "ov_pheno_GSE140082.qs"))

surv_data = select(
  pheno,
  "final_ostm:ch1",
  "final_osid:ch1"
) %>%
  rename("time" = 1, "status" = 2) %>%
  mutate_all(~ as.numeric(.))

# check
all(colnames(bulk) == rownames(surv_data))


set.seed(123)
alpha_samples <- Reduce(`*`, rep(c(2, 5), 2), init = 5e-4, accumulate = TRUE)

scab_obj <- create_scAB.v5(
  Object = seurat,
  bulk_dataset = bulk,
  phenotype = surv_data,
  method = 'survival',
  verbose = TRUE
)

k <- scAB::select_K.optimized(
  Object = scab_obj,
  K_max = 20L,
  repeat_times = 10L,
  maxiter = 2000L, # default in scAB
  seed = 123,
  verbose = TRUE
)

cli::cli_alert_success(crayon::green(
  "Optimal K = {k}"
))

scab_obj <- qs::qread(
  '/home/data/sigbridger/method_compare/survival/ov/GSE140082/scAB_obj.qs'
)
para_list <- select_alpha.optimized(
  Object = scab_obj,
  method = 'survival',
  K = k,
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
# 3093902
