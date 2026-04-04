library(data.table)
library(dplyr)
library(ScPP)
# library(BiocParallel)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/binary/brca/tnbc/GSE162228"
)

data_path = "/home/data/sigbridger/benchmark_data/brca"
# devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "seurat_tnbc.qs"))
search_res = fread("scpp_random_search.csv")

seurat_tum = readRDS(
  '/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_TNBCTum.rds'
)
tumor_cells = rownames(seurat_tum[[]])

search_res$benchmark = ifelse(colnames(seurat) %chin% tumor_cells, TRUE, FALSE)

# search_res = tibble::column_to_rownames(search_res, var = "V1")

source(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/compute_metrics.R"
)

metrics = compute_metrics(search_res)
cols <- rownames(metrics)
metrics = t(metrics)
colnames(metrics) = cols

# * random search, 51 times, from `scpp_search.R`
set.seed(162228)
probs_sample <- round(seq(0.01, 0.49, length.out = 49), 2)

arg_sample <- data.frame(
  prob = sample(probs_sample, 50, replace = TRUE),
  Log2FC_cutoff = round(runif(50), 3)
) %>%
  add_row(prob = 0.2, Log2FC_cutoff = 0.585) # default

probs = dplyr::bind_cols(arg_sample, metrics)

data.table::fwrite(probs, "scpp_acc.csv", row.names = T)
