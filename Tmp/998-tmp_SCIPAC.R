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

res <- DoSCIPAC(
  bulk,
  seurat,
  pheno,
  label_type = "SCIPAC",
  phenotype_class = "survival",
  hvg = 1000L,
  do_pca_sc = FALSE,
  n_pc = 60L,
  sc_batch_col = NULL,
  resolution = 2L,
  ela_net_alpha = 0.4,
  bt_size = 50L,
  ncore = 7L,
  ci_alpha = 0.05,
  nfold = 10L
)

qs::qsave(res, "/home/data/sigbridger/scipac_result.qs")

seurat <- seurat %>%
  SeuratObject::AddMetaData(rep("test", ncol(seurat)), col.name = "scissor") %>%
  AddMisc(
    scissor_type = "relapse",
    scissor_para = list(alpha = 0.05, cutoff = 0.2)
  )
