library(dplyr)
library(data.table)


# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/survival/brca/tnbc/TCGA_BRCA"
)

data_path = "/home/data/sigbridger/benchmark_data/brca"
devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "seurat_tnbc.qs"))

bulk = qs::qread(
  file.path(data_path, "brca_bulkdata_TCGA.qs")
)

pheno = qs::qread(file.path(
  data_path,
  "brca_surv_TCGA.qs"
))

surv_data = pheno

# bulk = bulk[, rownames(surv_data)]

if (!all(colnames(bulk) == rownames(surv_data))) {
  stop("bulk and surv_data not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}

# ! -------------------- scissor ---------------------

# * random search, 100 times
set.seed(123)


# * run scissor

scissor_result = Screen(
  matched_bulk = bulk,
  sc_data = seurat,
  phenotype = surv_data,
  label_type = "OS (M)_survival",
  phenotype_class = "survival",
  screen_method = "Scissor",
  alpha = NULL
)


pos_ratio = (scissor_result$scRNA_data$scissor == "Positive")
pos_ratio = as.data.frame(pos_ratio)
rownames(pos_ratio) = colnames(seurat)

seurat_tum = readRDS(
  '/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_TNBCTum.rds'
)
tumor_cells = rownames(seurat_tum[[]])
pos_ratio$benchmark = ifelse(colnames(seurat) %chin% tumor_cells, TRUE, FALSE)
source(
  '/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/compute_metrics.R'
)
res = compute_metrics(pos_ratio)

data.table::fwrite(
  res,
  file = "scissor_acc.csv",
  row.names = TRUE,
  col.names = TRUE
)
# 3060712
