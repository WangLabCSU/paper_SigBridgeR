library(data.table)
library(dplyr)
library(ScPP)
library(BiocParallel)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/survival/brca/her2/TCGA_BRCA"
)

data_path = "/home/data/sigbridger/benchmark_data/brca"

# * load data
seurat = qs::qread(file.path(data_path, "seurat_her2.qs"))

search_res = fread("scpp_random_search.csv")

seurat_tum = readRDS(
  '/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_HER2Tum.rds'
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
set.seed(8894)
probs = round(seq(0.01, 0.5, length.out = 49), 2)
probs = as.data.frame(probs)
probs = cbind(probs, metrics)

data.table::fwrite(probs, "scpp_acc.csv", row.names = T)
