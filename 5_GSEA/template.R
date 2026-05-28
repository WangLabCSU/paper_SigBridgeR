setwd(file.path(usethis::proj_path(), "5_GSEA"))

library(GSVA)


# * Hyperparameter
data_path <- "/home/data/sigbridger"
output_dir <- ""
param <- BiocParallel::MulticoreParam(workers = 4L)
seurat <- qs::qread(file.path(data_path, ""), nthreads = 8L)


methods <- c(
  "scissor",
  "scAB",
  "scPAS",
  "scPP",
  "SCIPAC",
  "DEGAS",
  "LP_SGL",
  "PIPET"
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

for (method in methods) {
  pos_cell <- which(seurat[[method]] == "Positive")
  non_pos_cell <- which(seurat[[method]] != "Positive")

  cli::cli_alert_info(
    "Positive cell: {length(pos_cell)}, non-Positive cell: {length(non_pos_cell)}"
  )

  expr_pos <- expr[, pos_cell]
  expr_non_pos <- expr[, non_pos_cell]

  gsea_param_pos <- gsvaParam(
    exprData = expr_pos,
    geneSets = gene_list,
    # kcdf = auto # * auto choose
  )
  gsea_param_non_pos <- gsvaParam(
    exprData = expr_non_pos,
    geneSets = gene_list,
    # kcdf = auto # * auto choose
  )

  gsea_res_pos <- gsva(
    gsea_param_pos,
    BPPARAM = param
  )
  gsea_res_non_pos <- gsva(
    gsea_param_non_pos,
    BPPARAM = param
  )

  data.table::fwrite(
    gsea_res_pos,
    file.path(output_dir, "gsea_pos.csv")
  )
  data.table::fwrite(
    gsea_res_non_pos,
    file.path(output_dir, "gsea_non_pos.csv")
  )

  cli::cli_h1("GSEA done")
}
