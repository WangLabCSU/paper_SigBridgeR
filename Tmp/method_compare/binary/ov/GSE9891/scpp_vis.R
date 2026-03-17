library(data.table)
library(dplyr)
library(ScPP)
library(BiocParallel)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/binary/ov/GSE9891"
)

data_path = "/home/data/sigbridger/benchmark_data/ov"
# devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "hgsoc_GSE165897_seurat.qs"))
search_res = fread("scpp_random_search.csv")

search_res$benchmark = ifelse(seurat$cell_type == 'EOC', TRUE, FALSE)

# search_res = tibble::column_to_rownames(search_res, var = "V1")

source(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/compute_metrics.R"
)

metrics = compute_metrics(search_res)
cols <- rownames(metrics)
metrics = t(metrics)
colnames(metrics) = cols
# * random search, 51 times, from `scpp_search.R`
set.seed(123)
probs_sample <- round(seq(0.01, 0.5, length.out = 49), 2)

arg_sample <- data.frame(
  prob = sample(probs_sample, 50, replace = TRUE),
  Log2FC_cutoff = round(runif(50), 3)
)

probs = dplyr::bind_cols(arg_sample, metrics)

data.table::fwrite(probs, "scpp_acc.csv", row.names = T)
