library(dplyr)
library(data.table)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/binary/brca/her2/GSE162228"
)

data_path = "/home/data/sigbridger/benchmark_data/brca"
devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "seurat_her2.qs"))

bulk = qs::qread(
  file.path(data_path, "brca_bulkdata_GSE162228.qs")
)

pheno = qs::qread(file.path(
  data_path,
  "brca_pheno_GSE162228.qs"
))

pheno_bi = setNames(
  case_when(
    pheno$`relapse status:ch1` == "relapse" ~ 1,
    pheno$`relapse status:ch1` == "non-relapse" ~ 0
  ),
  pheno$geo_accession
)


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
  label_type = "relapse",
  phenotype_class = "binary",
  screen_method = "Scissor",
  alpha = NULL,
  path2save_scissor_inputs = NULL
)


pos_ratio = (scissor_result$scRNA_data$scissor == "Positive")
pos_ratio = as.data.frame(pos_ratio)
rownames(pos_ratio) = colnames(seurat)

seurat_tum = readRDS(
  '/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_HER2Tum.rds'
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
# 720879
