# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  file.path(usethis::proj_path(), "4_positive_ctrl")
)
library(dplyr)
library(GSVA)
library(data.table)

data_path <- "brca"
markers_file_names <- "binary_deg_her2_GSE162228.csv"

# ? read marker file
raw_markers <- data.table::fread(file.path(data_path, markers_file_names))

# ? detect it contains at least 20 genes each
n_risk <- raw_markers[logFC > 0, .N]
n_prot <- raw_markers[logFC < 0, .N]
if (n_risk < 20 || n_prot < 20) {
  cli::cli_abort("Risk markers:  {n_risk},  Protective markers:  {n_prot}")
}


# ? convert to list and use top 20 (P.val < 0.05 first, log FC second)
top_risk <- raw_markers[logFC > 0][1:20] %>%
  dplyr::pull(gene)
top_protective <- raw_markers[logFC < 0][1:20] %>%
  dplyr::pull(gene)
gene_list <- list(
  "binary_GSE162228_pos_ssGSEA" = top_risk,
  "binary_GSE162228_neg_ssGSEA" = top_protective
)
# ? run ssGSEA
seurat_path <- "/home/data/sigbridger/benchmark_binary/brca/HER2"
seurat <- qs::qread(file.path(
  seurat_path,
  "binary_her2_GSE162228_merged_seurat.qs"
))

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
qs::qsave(
  es_df,
  file = "esmat/binary/brca/her2/GSE162228/ssGSEA_score_GSE162228.qs",
  nthreads = 4L
)
