# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/esmat/binary/ov/GSE140082"
)

library(GSVA)
library(dplyr)
library(data.table)

data_path <- "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/ov"
markers_file_names <- "binary_deg_GSE140082.csv"

# ? read marker file
raw_markers <- data.table::fread(file.path(data_path, markers_file_names))

# ? detect it contains at least 20 genes each
n_risk <- raw_markers[logFC > 0, .N]
n_prot <- raw_markers[logFC < 0, .N]
if (n_risk < 20 || n_prot < 20) {
  cli::cli_abort("Risk markers: ", n_risk, "; Protective markers: ", n_prot)
}

# ? convert to list and use top 20 (P.val < 0.05 first, log FC second)
top_risk <- raw_markers[logFC > 0][1:20] %>%
  dplyr::pull(gene)
top_protective <- raw_markers[logFC < 0][1:20] %>%
  dplyr::pull(gene)
gene_list <- list(
  "binary_GSE140082_pos_ssGSEA" = top_risk,
  "binary_GSE140082_neg_ssGSEA" = top_protective
)

# ? run ssGSEA
seurat_path <- "/home/data/sigbridger/benchmark_binary/ov/GSE165897"
seurat <- qs::qread(file.path(seurat_path, "GSE140082_ov_merged_seurat.qs"))

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
#   file = "ssGSEA_score_GSE140082.csv"
# )
qs::qsave(es_df, file = "ssGSEA_score_GSE140082.qs", nthreads = 4L)
# 43196

score = qs::qread("ssGSEA_score_GSE140082.qs", nthreads = 4L)
score = score[, 1:2]
score = cbind(score, seurat[[]])
qs::qsave(score, file = "ssGSEA_score_GSE140082.qs", nthreads = 4L)
