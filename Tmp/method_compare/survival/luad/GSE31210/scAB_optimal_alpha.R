library(dplyr)
library(Seurat)
library(BiocParallel)
library(scAB)
# library(SigBridgeR)

setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/survival/luad/GSE31210"
)

data_path = "/home/data/sigbridger/benchmark_data/lung"
devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "luad_GSE123902_seurat.qs"))

bulk = qs::qread(
  file.path(data_path, "lung_bulkdata_GSE31210.qs")
)

pheno = data.table::fread(file.path(
  "/home/data/data-resource/single-cell/Lung_Cancer/LUAD/lung.cancer.adeno.gse31210.hgu133plus2_entrezcdf.tsv"
))

surv_data = pheno %>%
  filter(`Array` %in% colnames(bulk)) %>%
  tibble::column_to_rownames("Array") %>%
  select("OS_Time", "OS_Status") %>%
  rename("time" = "OS_Time", "status" = "OS_Status")
bulk = bulk[, rownames(surv_data)]

if (!all(colnames(bulk) == rownames(surv_data))) {
  stop("bulk and surv_data not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}


# ! -------------------- scAB ---------------------

# * random search, 100 times
set.seed(123)
alpha_samples <- Reduce(`*`, rep(c(2, 5), 2), init = 5e-4, accumulate = TRUE)


scab_obj <- create_scAB.v5(
  Object = seurat,
  bulk_dataset = bulk,
  phenotype = surv_data,
  method = 'survival',
  verbose = TRUE
)

save_path = '/home/data/sigbridger/method_compare/survival/lung/GSE31210'
qs::qsave(scab_obj, file = file.path(save_path, "scAB_obj.qs"), nthreads = 4L)
# scab_obj <- qs::qread(file.path(save_path, 'scAB_obj.qs'))

k <- scAB::select_K.optimized(
  Object = scab_obj,
  K_max = 20L,
  repeat_times = 10L,
  maxiter = 2000L, # default in scAB
  seed = 123,
  verbose = TRUE
)

cli::cli_alert_success(crayon::green(
  "Optimal K: {k}"
))

para_list <- select_alpha.optimized(
  Object = scab_obj,
  method = 'survival',
  K = k,
  cross_k = 5,
  para_1_list = alpha_samples,
  para_2_list = alpha_samples,
  parallel = F
)

alpha <- para_list$para$alpha_1
alpha_2 <- para_list$para$alpha_2

cli::cli_alert_success(crayon::green(
  "Optimal parameters: alpha = {alpha}, alpha_2 = {alpha_2}"
))

qs::qsave(para_list, file = "scAB_para_list.qs", nthreads = 4L) # 3255601
