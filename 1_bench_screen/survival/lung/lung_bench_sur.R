# ! sc - GSE123902
# ! bulk - GSE3141, GSE8894, GSE31210, TCGA_LUAD
# ! survival
# ! lung (luad)

# ==============================================================================
# 1. Environment & Dependencies
# ==============================================================================
library(SigBridgeR)
library(Seurat)
library(dplyr)
library(rlang)
library(cli)
library(qs)

# 设置工作目录（建议后续改用 here:: 或 usethis::proj_path() 直接拼接绝对路径）
setwd(file.path(usethis::proj_path(), "1_bench_screen/survival/lung"))
data_path <- "/home/data/sigbridger/benchmark_data/lung"
save_path <- "/home/data/sigbridger/benchmark_survival/lung"

# 确保输出目录存在
dir.create(save_path, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. Configuration (集中管理数据集参数)
# ==============================================================================
bulk_configs <- list(
  GSE3141 = list(
    bulk_qs = "lung_bulkdata_GSE3141.qs",
    pheno_qs = "GSE3141_surv_pheno.qs",
    log_transform = TRUE,
    methods = c(
      "Scissor",
      "scAB",
      "SCIPAC",
      "scPAS",
      "scPP",
      "DEGAS",
      "LP_SGL",
      "PIPET"
    )
  ),
  GSE8894 = list(
    bulk_qs = "lung_bulkdata_GSE8894.qs",
    pheno_qs = "GSE8894_surv_pheno.qs",
    methods = c(
      "Scissor",
      "scAB",
      "SCIPAC",
      "scPAS",
      "scPP",
      "DEGAS",
      "LP_SGL",
      "PIPET"
    )
  ),
  GSE31210 = list(
    bulk_qs = "lung_bulkdata_GSE31210.qs",
    pheno_qs = "GSE31210_surv_pheno.qs",
    log_transform = TRUE,
    methods = c(
      "Scissor",
      "scAB",
      "SCIPAC",
      "scPAS",
      "scPP",
      "DEGAS",
      "LP_SGL",
      "PIPET"
    )
  ),
  TCGA_LUAD = list(
    bulk_qs = "TCGA_LUAD_bulkdata.qs",
    pheno_qs = "TCGA_LUAD_surv_pheno.qs",
    is_tcga = TRUE,
    log_transform = TRUE,
    methods = c(
      "Scissor",
      "scAB",
      "SCIPAC",
      "scPAS",
      "scPP",
      "DEGAS",
      "LP_SGL",
      "PIPET"
    )
  )
)

# ==============================================================================
# 3. Core Pipeline Function
# ==============================================================================
run_screening_pipeline <- function(
  config_name,
  config,
  sc_data,
  methods,
  data_path,
  save_path
) {
  cli::cli_h2("Starting pipeline for {.val {config_name}}")

  # 1. Load Bulk Data
  bulk <- qs::qread(file.path(data_path, config$bulk_qs), nthreads = 4)
  if (isTRUE(config$log_transform)) {
    bulk <- log2(bulk + 1)
  }

  # 2. Process Phenotype & Extract Survival Labels
  surv <- qs::qread(file.path(data_path, config$pheno_qs), nthreads = 4)

  # 3. Align Bulk Matrix with Labels
  cm_samples <- intersect(colnames(bulk), rownames(surv_data))
  bulk <- bulk[, cm_samples, drop = FALSE]
  surv_data <- surv_data[cm_samples, , drop = FALSE]
  cli::cli_alert_info(
    "Aligned bulk matrix: {nrow(bulk)} genes x {ncol(bulk)} samples"
  )

  # 4. Run Screening Methods
  results <- vector("list", length(methods))

  for (m in methods) {
    single_save_path <- file.path(
      save_path,
      paste0(config_name, "_", m, "_seurat.qs")
    )
    if (file.exists(single_save_path)) {
      cli::cli_alert_info(
        "Load result method: {.val {m}}, bulk: {.val {config_name}} (already exists)"
      )
      results[[m]] <- qs::qread(single_save_path, nthreads = 4L)
      next
    }

    screen_res <- SigBridgeR::Screen(
      bulk,
      sc_data,
      labels,
      label_type = paste0(m, "_binary"),
      phenotype_class = "binary",
      screen_method = m,
      alpha = if (m != "LP_SGL") NULL else 0.5,
      alpha_2 = NULL,
      path2save_scissor_inputs = NULL
    )

    # ? save directly for reproductivity and time saving
    qs::qsave(
      x = screen_res$scRNA_data, # a seurat object
      file = single_save_path
    )

    results[[m]] <- screen_res$scRNA_data
  }

  # 5. Merge & Save

  merged_res <- do.call(SigBridgeR::MergeResult, results)
  out_file <- file.path(
    save_path,
    paste0("survival_lung_", config_name, "_merged_seurat.qs")
  )
  qs::qsave(merged_res, out_file, nthreads = 8L)
  cli::cli_success("Saved to {.path {out_file}}\n")

  invisible(merged_res)
}

# ==============================================================================
# 4. Execution
# ==============================================================================
# 1. Load scRNA-seq data once
cli::cli_h1("Loading scRNA-seq reference...")
seurat_luad <- qs::qread(
  file.path(data_path, "luad_GSE123902_seurat.qs"),
  nthreads = 4
)

# 2. Set computational threads (TensorFlow & OpenMP)
SigBridgeR::setThreads(
  8L,
  tf_config = list(
    xla_flag = "--tf_xla_auto_jit=2 --tf_xla_cpu_global_jit",
    xla_device = NULL,
    inter_op = 8L,
    intra_op = 8L
  )
)

# 3. Run pipeline for all datasets sequentially
# (如需并行，可替换为 future.apply::future_lapply 或 parallel::mclapply)
lapply(names(bulk_configs), function(name) {
  if (
    file.exists(file.path(
      save_path,
      paste0("survival_lung_", name, "_merged_seurat.qs")
    ))
  ) {
    cli::cli_alert_info("Skipping {.val {name}} (already exists)")
    return(NULL)
  }
  run_screening_pipeline(
    name,
    bulk_configs[[name]],
    seurat_luad,
    bulk_configs[[name]]$methods,
    data_path,
    save_path
  )
})

cli::cli_h1("✅ All screening tasks completed.")
