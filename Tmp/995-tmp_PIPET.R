# setwd(dirname(.rs.api.getSourceEditorContext()$path))
setwd("~/R/Project/R_code/SigBridgeR/Tmp")
source("/home/yyx/R/Project/R_code/SigBridgeR/R/90-Test.R")

devtools::document("/home/yyx/R/Project/R_code/SigBridgeR")

data <- readRDS(
  "/home/yyx/R/Project/R_code/SigBridgeR/vignettes/example_data/binary_example_data.rds"
)

sc_mat <- data[[1]]
bulk <- data[[2]]
pheno <- data[[3]]

seurat <- SCPreProcess(
  sc_mat,
  meta_data = data.frame(
    test_col = rep("test", ncol(sc_mat)),
    row.names = colnames(sc_mat)
  ),
  # quality_control = FALSE,
  quality_control.pattern = c("^RP[LS]", "^MT-"),
  # data_filter = FALSE,
  data_filter.thresh = list(
    nFeature_RNA_thresh = c(200L, 6000L),
    # * only used when specifed in `quality_control.pattern`
    percent.mt = 20L, # mitochondrial genes
    percent.rp = 60L # ribosomal protein genes
    # ? When combined pattern is used, like `quality_control.pattern = "^MT-|^RP[LS]"`
    # ? Use `_` to separate different patterns like this:
    # percent.mt_rp = 60L

    # ? When filtering for non-mitochondrial genes and non-ribosomal proteins RNA genes,
    # ? the column names are in lowercase letter form with regular expression symbols removed.
    # `quality_control.pattern = "^[rt]rna"`
    # Correct threshhold setting is `percent.rt_rna = 60L`

    # ? Use `SigBridgeR:::Pattern2Colname` to get the correct colname if still confused.
  ),
  min_cells = 0,
  min_features = 0,
  scale_features = rownames(sc_mat),
  dims = 1:20,
  resolution = 0.1
)

pipet_result = Screen(
  matched_bulk = bulk,
  sc_data = seurat,
  phenotype = pheno,
  label_type = "binary_PIPET",
  phenotype_class = "binary",
  screen_method = "PIPET"
)

# devtools::document("/home/yyx/R/Project/R_code/SigBridgeR_methods/PIPET")

phenotype_df <- PIPET::AdaptPheno(
  phenotype = pheno,
  phenotype_type = "binary"
)

markers <- PIPET::Create_Markers2(
  bulk_data = bulk,
  colData = phenotype_df,
  class_col = "class"
)

pipet_result <- PIPET(
  sc_data = seurat,
  markers = markers,
  group = NULL,
  parallel = TRUE
)

data_path <- "/home/data/sigbridger"
qs::qsave(pipet_result, file.path(data_path, "pipet_result.qs"))

# 3362275
