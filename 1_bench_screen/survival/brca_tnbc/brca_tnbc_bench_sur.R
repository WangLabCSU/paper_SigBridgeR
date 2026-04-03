# ! sc - GSE161529
# ! bulk - GSE42568 , GSE162228 , TCGA_BRCA
# ! survival
# ! tnbc

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
setwd(file.path(usethis::proj_path(), "1_bench_screen/survival/brca_tnbc"))
data_path <- "/home/data/sigbridger/benchmark_data/brca"
save_path <- "/home/data/sigbridger/benchmark_survival/brca/TNBC"

# 确保输出目录存在
dir.create(save_path, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. Configuration (集中管理数据集参数)
# ==============================================================================
bulk_configs <- list(
  GSE42568 = list(
    bulk_qs = "brca_bulkdata_GSE42568.qs",
    pheno_qs = "brca_pheno_GSE42568.qs",
    id_col = "geo_accession",
    label_col = "tissue:ch1",
    label_map = c("breast cancer" = 1L, "normal breast" = 0L)
  ),
  GSE162228 = list(
    bulk_qs = "brca_bulkdata_GSE162228.qs",
    pheno_qs = "brca_pheno_GSE162228.qs",
    id_col = "geo_accession",
    label_col = "relapse status:ch1",
    label_map = c("relapse" = 1L, "non-relapse" = 0L)
  ),
  TCGA_BRCA = list(
    bulk_qs = "brca_bulkdata_TCGA.qs",
    pheno_qs = "brca_surv_TCGA.qs",
    is_tcga = TRUE,
    log_transform = TRUE
  )
)

methods <- c(
  "Scissor",
  "scAB",
  "SCIPAC",
  "scPAS",
  "scPP",
  "DEGAS",
  "LP_SGL",
  "PIPET"
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

  # 2. Process Phenotype & Extract Binary Labels
  pheno <- qs::qread(file.path(data_path, config$pheno_qs), nthreads = 4)

  if (config_name == "TCGA_BRCA") {
    surv_data <- pheno
  } else if (config_name == "GSE42568") {
    surv_data <- GSE42568_pheno %>%
      dplyr::select(
        `overall survival time_days:ch1`,
        `overall survival event:ch1`
      ) %>%
      dplyr::filter(
        !is.na(`overall survival time_days:ch1`) &
          !is.na(`overall survival event:ch1`) &
          `overall survival time_days:ch1` != "NA"
      ) %>%
      dplyr::rename(time = 1, status = 2)
  } else if (config_name == "GSE162228") {
    surv_data <- GSE162228_pheno %>%
      dplyr::select(`overall survival (years):ch1`, `alive:ch1`) %>%
      dplyr::rename(time = 1, status = 2) %>%
      dplyr::mutate(
        status = dplyr::case_when(
          status == "Alive" ~ 1L,
          status == "Death" ~ 0L
        )
      )
  }

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
    cli::cli_alert("Running {.val {m}}...")
    results[[m]] <- rlang::try_fetch(
      SigBridgeR::Screen(
        bulk,
        sc_data,
        surv_data,
        label_type = paste0(m, "_survival"),
        phenotype_class = "survival",
        screen_method = m,
        alpha = if (m != "LP_SGL") NULL else 0.5,
        alpha_2 = NULL,
        path2save_scissor_inputs = NULL
      ),
      error = function(e) {
        cli::cli_warn("{.fn {m}} failed: {.message {e$message}}")
        NULL
      }
    )
  }

  # 5. Merge & Save

  merged_res <- do.call(SigBridgeR::MergeResult, results)
  out_file <- file.path(
    save_path,
    paste0("survival_tnbc_", config_name, "_merged_seurat.qs")
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
seurat_tnbc <- qs::qread(file.path(data_path, "seurat_tnbc.qs"), nthreads = 4)

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
  run_screening_pipeline(
    name,
    bulk_configs[[name]],
    seurat_tnbc,
    methods,
    data_path,
    save_path
  )
})

cli::cli_h1("✅ All screening tasks completed.")
