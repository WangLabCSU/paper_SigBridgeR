# ! sc - GSE138852
# ! bulk - GSE28146 , GSE39420 , GSE109887
# ! binary
# ! ad

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
setwd(file.path(usethis::proj_path(), "1_bench_screen/binary/ad"))
data_path <- "/home/data/sigbridger/benchmark_data/ad"
save_path <- "/home/data/sigbridger/benchmark_binary/ad"

# 确保输出目录存在
dir.create(save_path, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. Configuration (集中管理数据集参数)
# ==============================================================================
bulk_configs <- list(
  GSE28146 = list(
    bulk_qs = "ad_bulkdata_GSE28146.qs",
    pheno_qs = "ad_pheno_GSE28146.qs",
    id_col = "geo_accession",
    label_col = "screen_label",
    methods = c(
      "Scissor",
      "scAB",
      #   "SCIPAC", # Caused by error in `sc.dat.rot[, 1:n.pc]`: subscript out of bounds
      "scPAS",
      "scPP",
      # "DEGAS", # IndexError: index 30 is out of bounds for axis 0 with size 30
      "LP_SGL"
      #   ,  "PIPET" # Warning: No markers found, please try different `lg2FC` and `p.adjust`
    )
  ),
  GSE39420 = list(
    bulk_qs = "ad_bulkdata_GSE39420.qs",
    pheno_qs = "ad_pheno_GSE39420.qs",
    id_col = "geo_accession",
    label_col = "screen_label",
    methods = c(
      "Scissor",
      "scAB",
      #   "SCIPAC", # sc.dat.rot[, 1:n.pc]
      "scPAS",
      "scPP",
      #   "DEGAS",
      "LP_SGL",
      "PIPET"
    )
  ),
  GSE109887 = list(
    bulk_qs = "ad_bulkdata_GSE109887.qs",
    pheno_qs = "ad_pheno_GSE109887.qs",
    id_col = "geo_accession",
    label_col = "screen_label",
    methods = c(
      "Scissor",
      "scAB",
      #   "SCIPAC",
      "scPAS",
      "scPP",
      #   "DEGAS",
      "LP_SGL"
      # ,      "PIPET" # Warning: ✖ Some classes have fewer than 2 marker genes, returning NULL
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

  # 2. Process Phenotype & Extract Binary Labels
  pheno <- qs::qread(file.path(data_path, config$pheno_qs), nthreads = 4)

  if (isTRUE(config$is_tcga)) {
    cm_samples <- intersect(pheno$sample, colnames(bulk))
    pheno <- pheno[pheno$sample %in% cm_samples, ]
    labels <- setNames(pheno$tumor, pheno$sample)
  } else if (!is.null(config$label_map)) {
    # 通用映射：利用命名向量索引，比 case_when 更高效且不易出错
    labels <- setNames(
      config$label_map[pheno[[config$label_col]]],
      pheno[[config$id_col]]
    )
    labels <- labels[!is.na(labels)] # 剔除未匹配的样本
  } else {
    labels <- setNames(pheno[[config$label_col]], pheno[[config$id_col]])
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
    if (m == "DEGAS") {
      cm_genes <- intersect(rownames(bulk), rownames(sc_data))
      bulk <- bulk[cm_genes, ]
      sc_data <- sc_data[cm_genes, ]
    } else if (file.exists(single_save_path)) {
      cli::cli_alert_info(
        "Load result method: {.val {m}}, bulk: {.val {config_name}} (already exists)"
      )
      results[[m]] <- qs::qread(single_save_path, nthreads = 4L)
      next
    }
    gc(verbose = FALSE)

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

  merged_res <- rlang::exec(SigBridgeR::MergeResult, !!!results)
  out_file <- file.path(
    save_path,
    paste0("binary_ad_", config_name, "_merged_seurat.qs")
  )
  qs::qsave(merged_res, out_file, nthreads = 8L)
  cli::cli_alert_success("Saved to {.path {out_file}}\n")

  invisible(merged_res)
}

# ==============================================================================
# 4. Execution
# ==============================================================================
# 1. Load scRNA-seq data once
cli::cli_h1("Loading scRNA-seq reference...")
seurat_ad <- qs::qread(
  file.path(data_path, "ad_GSE138852_seurat.qs"),
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

# 3. Run pipeline for all datasets
purrr::walk(names(bulk_configs), function(name) {
  if (
    file.exists(file.path(
      save_path,
      paste0("binary_ad_", name, "_merged_seurat.qs")
    )) &&
      !isTRUE(bulk_configs[[name]]$rerun)
  ) {
    cli::cli_alert_info("Skipping {.val {name}} (already exists)")
    return(NULL)
  }
  run_screening_pipeline(
    name,
    bulk_configs[[name]],
    seurat_ad,
    bulk_configs[[name]]$methods,
    data_path,
    save_path
  )
  NULL
})

cli::cli_h1("✅ All screening tasks completed.")
