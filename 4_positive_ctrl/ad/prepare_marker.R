setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(data.table)
library(Seurat)

# ! ad-sc-GSE138852
# ! 阳性对照组

data_path <- "/home/data/sigbridger/benchmark_data/ad"

b_bulks <- paste0("GSE", c("109887", "28146", "39420"))

# ! Single cell data
seurat <- qs::qread(
  file.path(data_path, "ad_GSE138852_seurat.qs"),
  nthreads = 4L
)
sc_genes <- rownames(seurat)

cli::cli_alert_info("Single cell data has {.val {length(sc_genes)}} genes")

cli::cli_h1("Finding marker genes for binary phenotype")

# ? Find markers for binary phenotype
purrr::walk(
  b_bulks,
  function(bulk_i) {
    bulk <- qs::qread(
      file.path(data_path, paste0("ad_bulkdata_", bulk_i, ".qs")),
      nthreads = 4L
    )
    bulk_genes <- rownames(bulk)
    cm_genes <- SeuratObject::intersect(sc_genes, bulk_genes)

    cli::cli_alert_info(
      "{bulk_i} has {.val {length(cm_genes)}} common genes with single cell data"
    )

    bulk <- bulk[cm_genes, ]

    pheno <- qs::qread(file.path(data_path, paste0("ad_pheno_", bulk_i, ".qs")))

    # 所有 pheno 均含现成的 screen_label (0/1) 列
    pheno <- pheno %>%
      dplyr::mutate(binary_screen = as.integer(screen_label)) %>%
      dplyr::select(geo_accession, binary_screen)

    if (any(bulk > 1000L)) {
      cli::cli_abort(
        "{bulk_i} has values greater than 1000, perhaps it is raw count matrix"
      )
    }

    if (!all(colnames(bulk) == rownames(pheno))) {
      cli::cli_alert_info("{bulk_i} matching pheno and bulk")
      cm_samples <- intersect(colnames(bulk), rownames(pheno))
      bulk <- bulk[, cm_samples]
      pheno <- pheno[cm_samples, ]
    }

    # Find markers
    group <- factor(
      pheno$binary_screen,
      levels = c(0, 1),
      labels = c("Control", "Case")
    )
    design <- stats::model.matrix(~group) # 截距模型：groupCase = Case - Control
    colnames(design) <- c("Intercept", "groupCase")

    fit <- limma::lmFit(bulk, design)
    fit <- limma::eBayes(fit, trend = TRUE) # 趋势化方差收缩，大样本推荐开启

    # 提取差异结果：Case vs Control
    deg <- limma::topTable(fit, coef = "groupCase", number = Inf, sort.by = "P")

    # 筛选显著 DEGs（常规阈值）
    sig_deg <- deg[abs(deg$logFC) >= 0.58 & deg$P.Value < 0.05, ]

    sig_deg <- sig_deg %>% dplyr::arrange(P.Value, dplyr::desc(abs(logFC)))

    if (nrow(sig_deg) < 50) {
      cli::cli_warn("{bulk_i} has less than 50 significant genes")
    }

    assign(paste0("binary_deg_", bulk_i), sig_deg, envir = .GlobalEnv)
  },
  .progress = "Binary pheno"
)

# ? Save marker genes of binary phenotype
purrr::walk(
  b_bulks,
  function(bulk_i) {
    data.table::fwrite(
      tibble::rownames_to_column(get(paste0("binary_deg_", bulk_i)), "gene"),
      file = paste0("binary_deg_ad_", bulk_i, ".csv")
    )
  }
)

cli::cli_h1("Done")
gc()
