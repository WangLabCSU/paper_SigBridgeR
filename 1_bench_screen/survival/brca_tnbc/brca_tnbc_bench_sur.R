# ! sc - GSE161529
# ! bulk - GSE42568, GSE162228, TCGA_BRCA
# ! survival
# ! TNBC

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
save_path <- "/home/data/sigbridger/benchmark_data/brca/TNBC"

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
    label_map = c("breast cancer" = 1L, "normal breast" = 0L),
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
  GSE162228 = list(
    bulk_qs = "brca_bulkdata_GSE162228.qs",
    pheno_qs = "brca_pheno_GSE162228.qs",
    id_col = "geo_accession",
    label_col = "relapse status:ch1",
    label_map = c("relapse" = 1L, "non-relapse" = 0L),
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
  TCGA_BRCA = list(
    bulk_qs = "brca_bulkdata_TCGA.qs",
    pheno_qs = "brca_pheno_TCGA.qs",
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

  # 2. Process Phenotype & Extract Binary Labels
  pheno <- qs::qread(file.path(data_path, config$pheno_qs), nthreads = 4)

  if (isTRUE(config$is_tcga)) {
    cm_samples <- intersect(pheno$sample, colnames(bulk))
    labels <- pheno %>%
      mutate(sample_type = substr(sample, 14, 15)) %>%
      filter(sample_type %in% c("01", "11"), sample %in% cm_samples) %>%
      mutate(sample_type = as.integer(sample_type == "01")) %>%
      {
        setNames(.$sample_type, .$sample)
      }
  } else {
    # 通用映射：利用命名向量索引，比 case_when 更高效且不易出错
    labels <- setNames(
      config$label_map[pheno[[config$label_col]]],
      pheno[[config$id_col]]
    )
    labels <- labels[!is.na(labels)] # 剔除未匹配的样本
  }

  # 3. Align Bulk Matrix with Labels
  bulk <- bulk[, names(labels), drop = FALSE]
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
      label_type = paste0(m, "_survival"),
      phenotype_class = "survival",
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
    paste0("survival_TNBC_", config_name, "_merged_seurat.qs")
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
seurat_TNBC <- qs::qread(file.path(data_path, "seurat_tnbc.qs"), nthreads = 4)

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
      paste0("survival_TNBC_", name, "_merged_seurat.qs")
    ))
  ) {
    cli::cli_alert_info("Skipping {.val {name}} (already exists)")
    return(NULL)
  }
  run_screening_pipeline(
    name,
    bulk_configs[[name]],
    seurat_TNBC,
    bulk_configs[[name]]$methods,
    data_path,
    save_path
  )
})

cli::cli_h1("✅ All screening tasks completed.")
