setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(data.table)
library(Seurat)
source("../../0_bench_data_preprocess/convert_gene_symbol.R")

# ! fshd-sc-GSE122873
# ! 阳性对照组

data_path <- "/home/data/sigbridger/benchmark_data/fshd"

b_bulks <- paste0("GSE", c("115650", "140261", "56787"))

# ? 各数据集原始 counts 文件（位于各自 GSE 文件夹中）
counts_files <- c(
  GSE115650 = "GSE115650/GSE115650_biopsy_genecounts.csv.gz",
  GSE140261 = "GSE140261/GSE140261_year2_genecounts.csv.gz",
  GSE56787 = "GSE56787/GSE56787_biopsy.gene.counts.csv.gz"
)

# ! Single cell data
seurat <- qs::qread(
  file.path(data_path, "fshd_GSE122873_seurat.qs"),
  nthreads = 4L
)
sc_genes <- rownames(seurat)

cli::cli_alert_info("Single cell data has {.val {length(sc_genes)}} genes")

cli::cli_h1("Finding marker genes for binary phenotype")

# ? Find markers for binary phenotype (raw counts + edgeR)
purrr::walk(
  b_bulks,
  function(bulk_i) {
    # 读取原始 counts，行名为 Ensembl ID
    counts <- data.table::fread(
      file.path(data_path, counts_files[[bulk_i]])
    )
    counts <- as.matrix(data.frame(counts, row.names = colnames(counts)[1]))
    storage.mode(counts) <- "numeric"

    # 去除 Ensembl 版本号，映射为 gene symbol 后与单细胞取交集
    ensembl <- sub("\\..*$", "", rownames(counts))
    sym <- convert_gene_symbol(ensembl)
    rownames(counts) <- sym

    # 同一 symbol 的多个 Ensembl 基因 counts 求和合并
    bulk <- rowsum(counts, group = sym)

    cm_genes <- intersect(sc_genes, rownames(bulk))

    cli::cli_alert_info(
      "{bulk_i} has {.val {length(cm_genes)}} common genes with single cell data"
    )

    bulk <- bulk[cm_genes, ]

    pheno <- qs::qread(
      file.path(data_path, paste0("fshd_pheno_", bulk_i, ".qs"))
    )

    if (bulk_i == "GSE115650") {
      col_names <- colnames(bulk) %>%
        stringr::str_remove("X") %>%
        stringr::str_replace("\\.", "-")
      sample_map <- setNames(pheno$geo_accession, pheno$title)
      colnames(bulk) <- sample_map[col_names]
    } else if (bulk_i == "GSE140261") {
      col_names <- colnames(bulk) %>%
        stringr::str_remove("X") %>%
        stringr::str_replace("\\.", "-")
      sample_map <- setNames(
        pheno$geo_accession,
        pheno$title %>% stringr::str_remove("_.*")
      )
      colnames(bulk) <- sample_map[col_names]
    } else if (bulk_i == "GSE56787") {
      col_names <- colnames(bulk)
      sample_map <- setNames(pheno$geo_accession, pheno$title)
      colnames(bulk) <- sample_map[col_names]
    }
    cm_samples <- intersect(colnames(bulk), rownames(pheno))
    bulk <- bulk[, cm_samples]
    pheno <- pheno[cm_samples, ]

    # Find markers with edgeR (quasi-likelihood pipeline)
    group <- factor(
      pheno$screen_label,
      levels = c(0, 1),
      labels = c("Control", "Case")
    )
    design <- model.matrix(~group) # 截距模型：groupCase = Case - Control
    colnames(design) <- c("Intercept", "groupCase")

    y <- edgeR::DGEList(counts = bulk, group = group)
    keep_expr <- edgeR::filterByExpr(y)
    y <- y[keep_expr, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y) # TMM 标准化
    y <- edgeR::estimateDisp(y, design)

    fit <- edgeR::glmQLFit(y, design, robust = TRUE)
    res <- edgeR::glmQLFTest(fit, coef = "groupCase")

    # 提取差异结果：Case vs Control
    deg <- edgeR::topTags(res, n = Inf, sort.by = "PValue")$table

    # 筛选显著 DEGs（常规阈值）
    sig_deg <- deg[abs(deg$logFC) >= 0.58 & deg$PValue < 0.05, ]

    sig_deg <- sig_deg %>% dplyr::arrange(PValue, desc(abs(logFC)))

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
      file = paste0("binary_deg_", bulk_i, ".csv")
    )
  }
)

gc()
