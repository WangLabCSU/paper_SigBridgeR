library(dplyr)
library(data.table)


# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/binary/brca/tnbc/GSE42568"
)

data_path = "/home/data/sigbridger/benchmark_data/brca"
devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "seurat_tnbc.qs"), nthreads = 4L)

bulk = qs::qread(
  file.path(data_path, "brca_bulkdata_GSE42568.qs")
)

pheno = qs::qread(file.path(
  data_path,
  "brca_pheno_GSE42568.qs"
))

pheno_bi = setNames(
  case_when(
    pheno$`tissue:ch1` == "breast cancer" ~ 1,
    pheno$`tissue:ch1` == "normal breast" ~ 0
  ),
  pheno$geo_accession
)

bulk = bulk[, names(pheno_bi)]

if (!all(colnames(bulk) == names(pheno_bi))) {
  stop("bulk and pheno_bi not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}

if (anyNA(pheno_bi)) {
  stop("pheno_bi has NA")
}


# ! -------------------- scissor ---------------------

# * random search, 100 times
set.seed(123)


# * run scissor

scissor_result = Screen(
  matched_bulk = bulk,
  sc_data = seurat,
  phenotype = pheno_bi,
  label_type = "tumor",
  phenotype_class = "binary",
  screen_method = "Scissor",
  alpha = NULL,
  path2save_scissor_inputs = NULL
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
# 1143547
