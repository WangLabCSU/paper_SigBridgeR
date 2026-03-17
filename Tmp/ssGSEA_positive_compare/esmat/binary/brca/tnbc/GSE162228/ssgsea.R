# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/esmat/binary/brca/tnbc/GSE162228"
)

library(GSVA)
library(dplyr)
library(data.table)

data_path <- "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/brca"
markers_file_names <- "binary_deg_tnbc_GSE162228.csv"

# ? read marker file
raw_markers <- data.table::fread(file.path(data_path, markers_file_names))

# ? detect it contains at least 20 genes each
n_risk <- raw_markers[logFC > 0, .N]
n_prot <- raw_markers[logFC < 0, .N]

# ! warning: the number of all markers are below 50
if (n_risk < 20 || n_prot < 20) {
  cli::cli_alert_warning(
    "Risk markers:  {n_risk},  Protective markers:  {n_prot}"
  )
}

# ? convert to list and use top 20 (P.val < 0.05 first, log FC second)
top_risk <- raw_markers[logFC > 0][1:20] %>%
  dplyr::pull(gene)
top_protective <- raw_markers[logFC < 0][1:20] %>%
  dplyr::pull(gene)
gene_list <- list(
  "binary_GSE162228_pos_ssGSEA" = top_risk[!is.na(top_risk)],
  "binary_GSE162228_neg_ssGSEA" = top_protective[!is.na(top_protective)]
)
# ? run ssGSEA
seurat_path <- "/home/data/sigbridger/benchmark_binary/brca/TNBC"
seurat <- qs::qread(file.path(seurat_path, "GSE162228_tnbc_merged_seurat.qs"))

expr <- as.matrix(SeuratObject::LayerData(
  seurat,
  layer = "data",
  assay = "RNA"
))

# ? run ssGSEA
param <- BiocParallel::MulticoreParam(workers = 2L)

ssgsea_param_sub <- gsvaParam(
  exprData = expr,
  geneSets = gene_list,
  # kcdf = auto # * auto choose
)
esmat_sub <- gsva(
  ssgsea_param_sub,
  BPPARAM = param
)

# ? convert to long form
es_df <- t(esmat_sub) %>% cbind(seurat[[]])

# ? save result
# data.table::fwrite(
#   es_df,
#   file = "ssGSEA_score_GSE162228.csv"
# )
qs::qsave(es_df, file = "ssGSEA_score_GSE162228.qs", nthreads = 4L)

score = qs::qread("ssGSEA_score_GSE162228.qs", nthreads = 4L)
score = score[, 1:2]
score = cbind(score, seurat[[]])
qs::qsave(score, file = "ssGSEA_score_GSE162228.qs", nthreads = 4L)
