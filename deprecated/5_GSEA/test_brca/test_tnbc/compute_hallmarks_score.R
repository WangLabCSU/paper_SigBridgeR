# ! TNBC

setwd(file.path(usethis::proj_path(), "5_GSEA"))

library(GSVA)


# * Hyperparameter
data_path <- "/home/data/sigbridger"
output_dir <- "/home/data/sigbridger/GSEA"
param <- BiocParallel::MulticoreParam(workers = 4L)

seurat_path <- "/home/data/sigbridger/benchmark_data/brca"
seurat <- qs::qread(
  file.path(
    seurat_path,
    "seurat_tnbc.qs"
  ),
  nthreads = 8L
)

geneset_hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
gene_list <- split(
  geneset_hallmark$gene_symbol,
  geneset_hallmark$gs_description
)

expr <- as.matrix(SeuratObject::LayerData(
  seurat,
  layer = "data",
  assay = "RNA"
))

gsea_param <- gsvaParam(
  exprData = expr,
  geneSets = gene_list,
  # kcdf = auto # * auto choose
)

gsea_res <- gsva(
  gsea_param,
  BPPARAM = param
)
data2save <- cbind(t(gsea_res), seurat[[]])

data.table::fwrite(data2save, file.path(output_dir, "tnbc_gsea_res.csv"))

cli::cli_h1("GSEA done")
