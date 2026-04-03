# ! sc - GSE165897
# ! bulk - GSE9891, GSE140082
# ! binary
# ! ov
# ! sc - GSE123902
# ! bulk - TCGA_LUAD
# ! binary
# ! LUAD

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
setwd(file.path(usethis::proj_path(), "1_bench_screen/binary/ov"))
data_path <- "/home/data/sigbridger/benchmark_data/ov"
save_path <- "/home/data/sigbridger/benchmark_binary/ov/"

# 确保输出目录存在
dir.create(save_path, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. Configuration (集中管理数据集参数)
# ==============================================================================
bulk_configs <- list(
  GSE9891 = list(
    bulk_qs = "ov_bulkdata_GSE9891.qs",
    pheno_qs = "ov_pheno_GSE9891.qs"
  ),
  GSE140082 = list(
    bulk_qs = "ov_bulkdata_GSE140082.qs",
    pheno_qs = "ov_pheno_GSE140082.qs"
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

  if (config_name == "GSE140082") {
    cli::cli_h1("GSE140082")
    pheno <- pheno[
      pheno$`newgrade:ch1` != "NA" &
        !is.na(pheno$`newgrade:ch1`),
    ]

    labels <- setNames(
      case_when(
        pheno$`newgrade:ch1` == "high.grade" ~ 1,
        pheno$`newgrade:ch1` == "low.grade" ~ 0
      ),
      pheno$geo_accession
    )

    # Handle NA in GSM4153781
    median_val <- median(bulk[, "GSM4153781"], na.rm = TRUE)
    na_indices <- which(is.na(bulk[, "GSM4153781"]))
    bulk[na_indices, "GSM4153781"] <- median_val
  } else if (config_name == "GSE9891") {
    cli::cli_h1("GSE9891")
    labels <- setNames(
      case_when(
        pheno$characteristics_ch1.1 == "Type : LMP" ~ 0,
        pheno$characteristics_ch1.1 == "Type : Malignant" ~ 1
      ),
      pheno$geo_accession
    )
  }

  # 3. Align Bulk Matrix with Labels
  bulk <- bulk[, names(labels), drop = FALSE]
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
        labels,
        label_type = paste0(m, "_binary"),
        phenotype_class = "binary",
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

  merged_res <- do.call(SigBridgeR::MergeResult, valid_results)
  out_file <- file.path(
    save_path,
    paste0("binary_ov_", config_name, "_merged_seurat.qs")
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
seurat_ov <- qs::qread(
  file.path(data_path, "hgsoc_GSE165897_seurat.qs"),
  nthreads = 8L
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
  run_screening_pipeline(
    name,
    bulk_configs[[name]],
    seurat_ov,
    methods,
    data_path,
    save_path
  )
})

cli::cli_h1("✅ All screening tasks completed.")
