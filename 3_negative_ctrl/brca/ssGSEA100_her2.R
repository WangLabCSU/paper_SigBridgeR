# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(file.path(usethis::proj_path(), "3_negative_ctrl/brca"))

set.seed(123)

library(GSVA)

# ! brca-her2

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_data/brca/seurat_her2.qs",
  nthreads = 4
)

expr <- as.matrix(SeuratObject::LayerData(seurat, layer = "data"))

# ! from Negative control markes
# ! random selected
markers <- data.table::fread(
  "her2_random20_markers_100rep.csv"
)

sample <- colnames(markers)

genesets <- as.list(markers)
names(genesets) <- paste0("Neg_ctrl_", sample)

ssgsea_param <- ssgseaParam(
  exprData = expr,
  geneSets = genesets
)

bp_param <- BiocParallel::MulticoreParam(workers = 4L)

mat_out <- gsva(ssgsea_param, BPPARAM = bp_param)

qs::qsave(
  mat_out,
  "her2_Sample_100_ssgsea_score.qs"
)
cli::cli_h1("Finish")
gc()
