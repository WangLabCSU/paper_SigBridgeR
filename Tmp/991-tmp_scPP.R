# ---- binary ----

library(ScPP)
load(system.file("data/binary.RData", package = "ScPP"))

sc = sc_Preprocess(sc_count)
geneList = marker_Binary(bulk, binary, ref_group = "Normal")
res = ScPP(sc, geneList)
Matrix::head(res$metadata)

# Phenotype+ genes
res$Genes_pos

# Phenotype- genes
res$Genes_neg

# Visualization of ScPP-identified cells
sc$ScPP = res$metadata$ScPP
SeuratObject::Idents(sc) = "ScPP"
DimPlot(sc, group = "ScPP", cols = c("grey", "blue", "red"))


# ---- continuous ----

library(ScPP)
load(system.file("data/continuous.RData", package = "ScPP"))

sc = sc_Preprocess(sc_count)
geneList = marker_Continuous(bulk, continuous$TMB_non_silent)
res = ScPP(sc, geneList)
Matrix::head(res$metadata)

# Phenotype+ genes
res$Genes_pos

# Phenotype- genes
res$Genes_neg

# Visualization of ScPP-identified cells
sc$ScPP = res$metadata$ScPP
SeuratObject::Idents(sc) = "ScPP"
DimPlot(sc, group = "ScPP", cols = c("grey", "blue", "red"))

# ---- survival ----

library(ScPP)
load(system.file("data/survival.RData", package = "ScPP"))

sc = ScPP::sc_Preprocess(sc_count)
geneList = marker_Survival(bulk, survival)
res = ScPP(sc, geneList)
Matrix::head(res$metadata)

# Phenotype+ genes
res$Genes_pos

# Phenotype- genes
res$Genes_neg

# Visualization of ScPP-identified cells
sc$ScPP = res$metadata$ScPP
SeuratObject::Idents(sc) = "ScPP"
DimPlot(sc, group = "ScPP", cols = c("grey", "blue", "red"))


# SigBridgeR 的示例数据
survival = select(survival, "time", "status")
write.csv(
  bulk,
  "/home/yyx/R/Project/R_code/SigBridgeR/vignettes/example_data/survival_bulk.csv"
)
write.csv(
  survival,
  "/home/yyx/R/Project/R_code/SigBridgeR/vignettes/example_data/survival_survival.csv"
)

saveRDS(
  object = list(mat_exam, bulk, survival),
  file = "/home/yyx/R/Project/R_code/SigBridgeR/vignettes/example_data/survival_example_data.rds"
)

# 复现
sc = SCPreProcess(
  sc_count
)
geneList = marker_Survival(bulk, survival)
geneList = marker_Survival(bulk, pheno)

res = ScPP(seurat, geneList)

seurat = ScPP::sc_Preprocess(mat_exam)


scpp_result = Screen(
  bulk,
  seurat,
  survival,
  label_type = "survival_scPP",
  phenotype_class = "survival",
  screen_method = "scPP"
)
scpp_result = Screen(
  bulk,
  sc,
  survival,
  label_type = "survival_scPP",
  phenotype_class = "survival",
  screen_method = "scPP"
)


scpp_result = Screen(
  bulk,
  seurat,
  pheno,
  label_type = "survival_scPP",
  phenotype_class = "survival",
  screen_method = "scPP",
  probs = NULL
)
