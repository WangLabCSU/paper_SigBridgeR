# setwd(dirname(.rs.api.getSourceEditorContext()$path))
setwd("~/R/Project/R_code/SigBridgeR/Tmp")

devtools::document("/home/yyx/R/Project/R_code/SigBridgeR")

data <- readRDS(
  "/home/yyx/R/Project/R_code/SigBridgeR/vignettes/example_data/survival_example_data.rds"
)

bulk <- data[[2]]
pheno <- data[[3]]

seurat <- qs::qread("/home/data/sigbridger/merged_without_pipet.qs")
# seurat <- seurat$scRNA_data

res <- DoSIDISH(bulk, seurat, pheno, sidish_param = list(device = 'cpu'))

qs::qsave(res, "/home/data/sigbridger/sidish_result.qs")
