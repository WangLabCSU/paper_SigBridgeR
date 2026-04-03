setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd(
#   "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/ov"
# )

library(dplyr)
library(data.table)

# ! ov-sc-GSE165897
# ! 阳性对照组

data_path <- "/home/data/sigbridger/benchmark_data/lung"

b_bulks <- c("TCGA_LUAD")
s_bulks <- c("TCGA_LUAD", paste0("GSE", c("3141", "8894", "31210")))

# ! Single cell data
seurat <- qs::qread(
  file.path(data_path, "luad_GSE123902_seurat.qs"),
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
      file.path(data_path, paste0(bulk_i, "_bulkdata.qs")),
      nthreads = 4L
    )

    bulk_genes <- rownames(bulk)
    cm_genes <- intersect(sc_genes, bulk_genes)

    cli::cli_alert_info(
      "{bulk_i} has {.val {length(cm_genes)}} common genes with single cell data"
    )

    bulk <- bulk[cm_genes, ]

    pheno <- qs::qread(file.path(
      data_path,
      paste0(bulk_i, "_pheno.qs")
    ))

    if (bulk_i == "TCGA_LUAD") {
      pheno <- mutate(pheno, sample_type = substr(pheno$sample, 14, 15)) %>%
        select(sample, sample_type) %>%
        filter(sample_type %in% c("01", "11")) %>%
        mutate(sample_type = ifelse(sample_type == "01", 1, 0))
    }

    if (nrow(pheno) != ncol(bulk)) {
      cli::cli_alert_info("{bulk_i} matching pheno and bulk")
      cm_samples <- intersect(colnames(bulk), pheno$sample)
      bulk <- bulk[, cm_samples]
      pheno <- pheno[pheno$sample %in% cm_samples, ]
    }

    if (any(bulk > 1000)) {
      cli::cli_warn(
        "{bulk_i} has values greater than 1000, perhaps it is raw count matrix"
      )

      dge <- edgeR::DGEList(counts = bulk, group = pheno$sample_type) # counts: integer matrix
      dge <- edgeR::calcNormFactors(dge) # TMM normalization

      # 2. voom 转换：counts → log-CPM + weights
      design <- model.matrix(
        ~ factor(
          pheno$sample_type,
          levels = c(0, 1),
          labels = c("Control", "Case")
        )
      )
      colnames(design) <- c("Intercept", "groupCase")

      v <- limma::voom(
        dge,
        design = design,
        plot = TRUE
      ) # design = model.matrix(~ group)

      fit <- limma::lmFit(v, design)
      fit <- limma::eBayes(fit)
    } else {
      # Find markers
      group <- factor(
        pheno$binary_screen,
        levels = c(0, 1),
        labels = c("Control", "Case")
      )
      design <- model.matrix(~group) # 截距模型：groupCase = Case - Control
      colnames(design) <- c("Intercept", "groupCase")

      fit <- limma::lmFit(bulk, design)
      fit <- limma::eBayes(fit, trend = TRUE) # 趋势化方差收缩，大样本推荐开启
    }

    # 4. 提取差异结果：Case vs Control
    deg <- limma::topTable(fit, coef = "groupCase", number = Inf, sort.by = "P")

    # 5. 筛选显著 DEGs（常规阈值）
    sig_deg <- deg[abs(deg$logFC) >= 0.58 & deg$P.Value < 0.05, ]

    sig_deg <- sig_deg %>% dplyr::arrange(P.Value, desc(abs(logFC)))

    if (nrow(sig_deg) < 50) {
      cli::cli_warn("{bulk_i} has less than 50 significant genes")
    }

    # library(ggrepel)
    # deg$Sig <- ifelse(deg$adj.P.Val < 0.05 & abs(deg$logFC) >= 1, "Yes", "No")
    # ggplot(deg, aes(logFC, -log10(P.Value), color = Sig)) +
    #   geom_point(alpha = 0.6) +
    #   scale_color_manual(values = c("No" = "gray30", "Yes" = "red")) +
    #   theme_minimal() +
    #   xlab("log2 Fold Change") +
    #   ylab("-log10(P-value)")

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

# ? Find markers for survival phenotype
future::plan(future::multicore, workers = 8)
cli::cli_h1("Finding marker genes for survival phenotype")

purrr::walk(
  s_bulks,
  function(bulk_i) {
    bulk <- if (bulk_i == "TCGA_LUAD") {
      qs::qread(file.path(data_path, "TCGA_LUAD_bulkdata.qs"), nthreads = 4L)
    } else {
      qs::qread(
        file.path(data_path, paste0("lung_bulkdata_", bulk_i, ".qs")),
        nthreads = 4L
      )
    }

    bulk_genes <- rownames(bulk)
    cm_genes <- intersect(sc_genes, bulk_genes)

    cli::cli_alert_info(
      "{bulk_i} has {.val {length(cm_genes)}} common genes with single cell data"
    )

    bulk <- bulk[cm_genes, ]

    if (bulk_i %chin% c("GSE3141", "GSE8894", "GSE31210")) {
      surv_path <- if (bulk_i == "GSE3141") {
        "lung.cancer.adeno.gse3141.hgu133plus2_entrezcdf.tsv"
      } else if (bulk_i == "GSE8894") {
        "lung.cancer.adeno.gse8894.gpl570.tsv"
      } else if (bulk_i == "GSE31210") {
        "lung.cancer.adeno.gse31210.hgu133plus2_entrezcdf.tsv"
      }

      surv <- data.table::fread(file.path(
        "/home/data/data-resource/single-cell/Lung_Cancer/LUAD/",
        surv_path
      )) %>%
        tibble::column_to_rownames("Array") %>%
        dplyr::select("OS_Time", "OS_Status") %>%
        dplyr::rename("time" = "OS_Time", "status" = "OS_Status")
    }

    if (bulk_i == "TCGA_LUAD") {
      surv <- qs::qread(file.path(data_path, "TCGA_LUAD_surv_pheno.qs"))
    }

    if (any(rownames(surv) != colnames(bulk))) {
      cli::cli_alert_info("In {bulk_i}: matching pheno and surv")
      cm_samples <- intersect(colnames(bulk), rownames(surv))
      bulk <- bulk[, cm_samples]
      surv <- surv[cm_samples, ]
    }

    # remove duplicated genes
    dt <- data.table::as.data.table(bulk, keep.rownames = "gene")
    dt <- dt[, lapply(.SD, max), by = gene, .SDcols = names(dt)[-1]]
    bulk <- as.data.frame(dt) %>% tibble::column_to_rownames("gene") %>% t()

    if (!all(rownames(surv) == rownames(bulk))) {
      cli::cli_abort("In {bulk_i}: final matching failed")
    }

    if (any(bulk > 1000)) {
      cli::cli_warn(
        "{bulk_i} has values greater than 1000, perhaps it is raw count matrix, log2 transforming"
      )
      bulk <- log2(bulk + 1L)
    }

    # cox
    surv_obj <- survival::Surv(surv$time, surv$status)

    cox_results <- furrr::future_map_dfr(
      colnames(bulk),
      ~ {
        x <- as.numeric(bulk[, .x])
        fit <- survival::coxph(surv_obj ~ x)
        s <- Matrix::summary(fit)
        tibble::tibble(
          gene = .x,
          beta = stats::coef(fit),
          se = sqrt(stats::vcov(fit)[1, 1]),
          pval = s$sctest[3], # LRT p-value
          HR = exp(stats::coef(fit)),
          HR_low = exp(stats::coef(fit) - 1.96 * sqrt(stats::vcov(fit)[1, 1])),
          HR_up = exp(stats::coef(fit) + 1.96 * sqrt(stats::vcov(fit)[1, 1]))
        )
      },
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )

    cox_results <- cox_results %>%
      mutate(
        FDR = p.adjust(pval, method = "BH"),
        sig = (pval < 0.05) & (!is.na(pval)),
        direction = case_when(
          HR > 1 ~ "risk",
          HR < 1 ~ "protective",
          TRUE ~ "background"
        ),
        abs_logHR = abs(log(HR))
      ) %>%
      arrange(pval, desc(abs_logHR))

    filtered_cox_results <- cox_results %>% # 按 HR 降序 → 顶部即「高 HR 基因」
      dplyr::filter(sig, !is.na(HR))

    if (nrow(filtered_cox_results) < 50) {
      cli::cli_h2("{bulk_i} has less than 50 significant genes in survival")
    } else if (nrow(filtered_cox_results) == 0) {
      cli::cli_h2(cli::col_red("{bulk_i} has no significant genes in survival"))
      return(NULL)
    }

    assign(
      paste0("surv_cox_", bulk_i),
      filtered_cox_results,
      envir = .GlobalEnv
    )
  },
  .progress = "Survival pheno"
)

# ? Save marker genes of survival phenotype
purrr::walk(
  s_bulks,
  function(bulk_i) {
    # save
    data.table::fwrite(
      get(paste0("surv_cox_", bulk_i)),
      file = paste0("survival_deg_", bulk_i, ".csv")
    )
  }
)
