# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/esmat/survival/luad/GSE31210"
)

library(dplyr)
library(GSVA)
library(data.table)

data_path <- "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/luad"
markers_file_names <- "survival_deg_GSE31210.csv"

# ? read marker file
raw_markers <- data.table::fread(file.path(data_path, markers_file_names))

# ? detect it contains at least 20 genes each
n_risk <- raw_markers[direction == "risk", .N]
n_prot <- raw_markers[direction == "protective", .N]
if (n_risk < 20 || n_prot < 20) {
  cli::cli_abort("Risk markers: ", n_risk, "; Protective markers: ", n_prot)
}

# ? convert to list and use top 20
top_risk <- raw_markers[direction == "risk"][order(-abs_logHR)][1:20] %>%
  dplyr::pull(gene)
top_protective <- raw_markers[direction == "protective"][order(-abs_logHR)][
  1:20
] %>%
  dplyr::pull(gene)
gene_list <- list(
  "survival_GSE31210_pos_ssGSEA" = top_risk,
  "survival_GSE31210_neg_ssGSEA" = top_protective
)

# ? run ssGSEA
seurat_path <- "/home/data/sigbridger/benchmark_data/lung/GSE31210"
seurat <- qs::qread(file.path(seurat_path, "GSE31210_luad_merged_seurat.qs"))

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
#   file = "ssGSEA_score_GSE31210.csv"
# )
qs::qsave(es_df, file = "ssGSEA_score_GSE31210.qs", nthreads = 4L)
# 3995965

# score = qs::qread("ssGSEA_score_GSE31210.qs")
# diff_cols <- grepv("sc|DEGAS|LP_SGL|PIPET", colnames(seurat[[]]))
# score2 <- cbind(score[, c(1, 2)], seurat[[]])
# qs::qsave(score2, file = "ssGSEA_score_GSE31210.qs", nthreads = 4L)

# raw <- qs::qread(
#   "/home/data/sigbridger/benchmark_data/lung/luad_GSE123902_seurat.qs"
# )
