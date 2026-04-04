library(dplyr)
library(data.table)


# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/binary/ov/GSE9891"
)

data_path = "/home/data/sigbridger/benchmark_data/ov"
devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "hgsoc_GSE165897_seurat.qs"))

bulk = qs::qread(
  file.path(data_path, "ov_bulkdata_GSE9891.qs"),
  nthreads = 4
)

pheno = qs::qread(file.path(data_path, "ov_pheno_GSE9891.qs"))

pheno_bi = setNames(
  case_when(
    pheno$characteristics_ch1.1 == "Type : LMP" ~ 0,
    pheno$characteristics_ch1.1 == "Type : Malignant" ~ 1
  ),
  pheno$geo_accession
)

# ! 匹配
bulk = bulk[, names(pheno_bi)]

if (any(colnames(bulk) != names(pheno_bi))) {
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
  label_type = "grade",
  phenotype_class = "binary",
  screen_method = "Scissor",
  alpha = NULL,
  path2save_scissor_cache = NULL
)


pos_ratio = (scissor_result$scRNA_data$scissor == "Positive")
pos_ratio = as.data.frame(pos_ratio)
rownames(pos_ratio) = colnames(seurat)

pos_ratio$benchmark = ifelse(seurat$cell_type == 'EOC', TRUE, FALSE)
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
# 3435881
# 381754
