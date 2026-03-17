# setwd(dirname(.rs.api.getActiveDocumentContext()$path))
setwd("~/R/Project/R_code/SigBridgeR/Tmp/GSEA/ov")

set.seed(123)

library(Seurat)
# library(clusterProfiler)
# library(dplyr)
library(ggplot2)
library(RcppML)
library(irGSEA)
# devtools::document("/home/yyx/R/Project/R_code/SigBridgeR")

# ! bulk-GSE140082, sc-GSE165897, survival as phenotype
seurat <- qs::qread(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_data/ov/GSE165897/GSE140082_ov_merged_seurat.qs",
  nthreads = 4
)

# SaveSeuratRds(seurat, file = "GSE140082_ov_merged_seurat.rds")

# seurat = SCPreProcess(seurat, verbose = FALSE) # * update seurat

# ! 速度慢
seurat_score <- irGSEA.score(
  object = seurat,
  assay = "RNA",
  slot = "data",
  seeds = 123,
  ncores = 4,
  min.cells = 3,
  min.feature = 0,
  custom = F,
  geneset = NULL,
  msigdb = T,
  species = "Homo sapiens",
  category = "C5",
  subcategory = NULL,
  geneid = "symbol",
  method = c("AUCell", "UCell", "singscore", "ssgsea", "JASMINE", "viper"),
  aucell.MaxRank = NULL,
  ucell.MaxRank = NULL,
  kcdf = 'Gaussian'
)

qs::qsave(
  seurat_score,
  "/home/data/sigbridger/GSEA/ov/ov_GSE140082_irGSEA_score.qs",
  nthreads = 4L
)
# PID=89334
