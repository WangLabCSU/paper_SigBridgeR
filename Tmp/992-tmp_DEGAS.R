# setwd(dirname(.rs.api.getSourceEditorContext()$path))
setwd("~/R/Project/R_code/SigBridgeR/Tmp")
# source("/home/yyx/R/Project/R_code/SigBridgeR/R/90-Test.R")

devtools::document("/home/yyx/R/Project/R_code/SigBridgeR")
save_path = '/home/data/sigbridger'

data <- readRDS(
  "/home/yyx/R/Project/R_code/SigBridgeR/vignettes/example_data/survival_example_data.rds"
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

seurat <- qs::qread(
  "/home/data/sigbridger/merged_without_pipet.qs",
  nthreads = 4L
)

# ListPyEnv()
# Warning: No venv found in "~/.virtualenvs", "~/.venvs", "./venv", and "./.venv", return empty result
#     name                                      python  type
# 1   base             /home/yyx/miniconda3/bin/python conda
# 2 degas2 /home/yyx/miniconda3/envs/degas2/bin/python conda

# SetupPyEnv(
#     env_type = "conda",
#     env_name = "r-reticulate-degas",
#     method = c("system"),
#     env_file = "/home/yyx/R/Project/R_code/SigBridgeR/inst/conda/DEGAS_environment.yml",
#     recreate = TRUE,
#     use_conda_forge = TRUE,
#     verbose = TRUE
# )

# # ! 不必要的，DODEGAS里集成了
# reticulate::use_condaenv(
#     "/home/yyx/miniconda3/envs/r-reticulate-degas/bin/python"
# )

# debug(runCCMTL.optimized)

# degas_result <- Screen(
#   bulk,
#   seurat,
#   pheno,
#   label_type = "This_is_a_DEGAS_test",
#   phenotype_class = "survival",
#   screen_method = "DEGAS"
# )

# qs::qsave(degas_result, "/home/data/sigbridger/degas_result.qs") # 3450855

scissor_result = Screen(
  path2load_scissor_cache = "/home/yyx/R/Project/R_code/SigBridgeR/Scissor_inputs.RData",
  matched_bulk = bulk,
  sc_data = seurat,
  phenotype = pheno,
  label_type = "survival",
  phenotype_class = "survival",
  screen_method = "Scissor",
  alpha = 0.05,
  reliability_test = list(
    run = TRUE,
    n = 10L,
    nfold = 10L
  )
)
qs::qsave(
  scissor_result,
  file.path(save_path, "scissor_result.qs"),
  nthreads = 4
)

# scpas_result = Screen(
#     matched_bulk = bulk,
#     sc_data = seurat,
#     phenotype = pheno,
#     label_type = "survival",
#     phenotype_class = "survival",
#     screen_method = "scPAS",
#     alpha = NULL
# )
# qs::qsave(scpas_result, file.path(save_path, "scpas_result.qs"), nthreads = 4)

# scab_result = Screen(
#     matched_bulk = bulk,
#     sc_data = seurat,
#     phenotype = pheno,
#     label_type = "survival",
#     phenotype_class = "survival",
#     screen_method = "scAB",
#     alpha = 0.005,
#     alpha_2 = 0.005
# )
# qs::qsave(scab_result, file.path(save_path, "scab_result.qs"), nthreads = 4)

# scpp_result = Screen(
#     matched_bulk = bulk,
#     sc_data = seurat,
#     phenotype = pheno,
#     label_type = "survival",
#     phenotype_class = "survival",
#     screen_method = "scPP",
#     probs = 0.2
# )

# qs::qsave(scpp_result, file.path(save_path, "scpp_result.qs"), nthreads = 4) # 3594499

# SaveEnv()

# # ScreenFractionPlot(scissor_result$scRNA_data, group_by = "seurat_clusters")

# LoadEnv()

# ScreenFractionPlot(scissor_result$scRNA_data, group_by = "seurat_clusters")
# merged_seurat <- MergeResult(
#     scissor_result$scRNA_data,
#     scpas_result$scRNA_data,
#     scab_result$scRNA_data,
#     scpp_result$scRNA_data
# )
# fraction = ScreenFractionPlot(merged_seurat, group_by = "seurat_clusters")
# upset = ScreenUpsetPlot(merged_seurat)
