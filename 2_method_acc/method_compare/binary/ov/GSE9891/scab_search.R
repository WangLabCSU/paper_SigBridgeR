library(dplyr)
library(Seurat)
# library(BiocParallel)
library(scAB)
# library(SigBridgeR)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/binary/ov/GSE9891"
)
data_path = "/home/data/sigbridger/benchmark_data/ov"

# * load data
seurat = qs::qread(file.path(data_path, "hgsoc_GSE165897_seurat.qs"))

# ! -------------------- scAB ---------------------

# * random search, 100 times
set.seed(123)
alpha_samples <- Reduce(`*`, rep(c(2, 5), 2), init = 5e-4, accumulate = TRUE)


save_path = '/home/data/sigbridger/method_compare/binary/ov/GSE9891'
# qs::qsave(scab_obj, file = file.path(save_path, "scAB_obj.qs"), nthreads = 4L)
scab_obj <- qs::qread(file.path(save_path, "scab_obj.qs"))


scab_res = scAB::scAB.optimized(
  Object = scab_obj,
  K = 2L,
  alpha = 0.05,
  alpha_2 = 0.005
)

seurat_screened = scAB::findSubset.optimized(
  Object = seurat,
  scAB_Object = scab_res,
  tred = 2L
)

# 合并所有结果
all_results <- ifelse(seurat_screened$scAB == 'Positive', TRUE, FALSE)
all_results <- as.data.frame(all_results)
rownames(all_results) = colnames(seurat) # each cell is a row

all_results$benchmark = ifelse(seurat$cell_type == 'EOC', TRUE, FALSE)

data.table::fwrite(
  all_results,
  file = "scab_random_search.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("scAB random search completed."))

# ! ----- Visualization ----

# * 计算指标
source(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/compute_metrics.R"
)
metrics = compute_metrics(all_results)

data.table::fwrite(
  metrics,
  file = "scab_acc.csv",
  row.names = TRUE
)

# # PID=3099497
