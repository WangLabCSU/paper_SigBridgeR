setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd(
#   "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/ov"
# )

library(dplyr)
library(data.table)

# ! ov-sc-GSE165897
# ! 阳性对照组

data_path <- "/home/data/sigbridger/benchmark_data/brca"

b_bulks <- c("TCGA", paste0("GSE", c("42568", "162228")))
s_bulks <- c("TCGA", paste0("GSE", c("42568", "162228")))

# ! Single cell data
seurat_her2 <- qs::qread(
  file.path(data_path, "seurat_her2.qs"),
  nthreads = 4L
)
sc_genes_her2 <- rownames(seurat_her2)

cli::cli_alert_info(
  "HER2: Single cell data has {.val {length(sc_genes_her2)}} genes"
)

seurat_tnbc <- qs::qread(
  file.path(data_path, "seurat_tnbc.qs"),
  nthreads = 4L
)
sc_genes_tnbc <- rownames(seurat_tnbc)
cli::cli_alert_info(
  "TNBC: Single cell data has {.val {length(sc_genes_her2)}} genes"
)


purrr::walk(c("her2", "tnbc"), function(subtype) {
  sc_genes <- get(paste0("sc_genes_", subtype))

  cli::cli_h1("Processing subtype: {subtype}")

  cli::cli_h2("Finding marker genes for binary phenotype")

  # ? Find markers for binary phenotype
  purrr::walk(
    b_bulks,
    function(bulk_i) {
      bulk <- qs::qread(
        file.path(data_path, paste0("brca_bulkdata_", bulk_i, ".qs")),
        nthreads = 4L
      )
      rn <- rownames(bulk)
      rn[is.na(rn)] <- paste0("unknown_", which(is.na(rn)))
      rownames(bulk) <- rn

      bulk_genes <- rownames(bulk)
      cm_genes <- intersect(sc_genes, bulk_genes)

      cli::cli_alert_info(
        "{bulk_i} has {.val {length(cm_genes)}} common genes with single cell data"
      )

      bulk <- bulk[cm_genes, ]

      pheno <- qs::qread(file.path(
        data_path,
        paste0("brca_pheno_", bulk_i, ".qs")
      ))

      if (bulk_i == "GSE42568") {
        pheno <- pheno %>%
          {
            setNames(
              case_when(
                .$`tissue:ch1` == "breast cancer" ~ 1,
                .$`tissue:ch1` == "normal breast" ~ 0
              ),
              .$geo_accession
            )
          }
      } else if (bulk_i == "GSE162228") {
        pheno <- pheno %>%
          {
            setNames(
              case_when(
                .$`relapse status:ch1` == "relapse" ~ 1,
                .$`relapse status:ch1` == "non-relapse" ~ 0
              ),
              .$geo_accession
            )
          }
      } else if (bulk_i == "TCGA") {
        pheno <- dplyr::mutate(
          pheno,
          sample_type = substr(pheno$sample, 14, 15)
        ) %>%
          dplyr::select(sample, sample_type) %>%
          dplyr::filter(
            sample_type %in% c("01", "11")
          ) %>%
          dplyr::mutate(sample_type = ifelse(sample_type == "01", 1, 0))

        pheno <- setNames(pheno$sample_type, pheno$sample)
      }

      if (length(pheno) != ncol(bulk)) {
        cli::cli_alert_info("{bulk_i} matching pheno and bulk")
        cm_samples <- intersect(colnames(bulk), names(pheno))
        bulk <- bulk[, cm_samples]
        pheno <- pheno[names(pheno) %in% cm_samples]
      }

      if (any(bulk > 1000)) {
        cli::cli_warn(
          "{bulk_i} has values greater than 1000, perhaps it is raw count matrix. using voom"
        )

        dge <- edgeR::DGEList(counts = bulk, group = pheno) # counts: integer matrix
        dge <- edgeR::calcNormFactors(dge) # TMM normalization

        # 2. voom 转换：counts → log-CPM + weights
        design <- model.matrix(
          ~ factor(
            pheno,
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
          pheno,
          levels = c(0, 1),
          labels = c("Control", "Case")
        )
        design <- model.matrix(~group) # 截距模型：groupCase = Case - Control
        colnames(design) <- c("Intercept", "groupCase")

        fit <- limma::lmFit(bulk, design)
        fit <- limma::eBayes(fit, trend = TRUE) # 趋势化方差收缩，大样本推荐开启
      }

      # 4. 提取差异结果：Case vs Control
      deg <- limma::topTable(
        fit,
        coef = "groupCase",
        number = Inf,
        sort.by = "P"
      )

      # 5. 筛选显著 DEGs（常规阈值）
      sig_deg <- deg[abs(deg$logFC) >= 0.58 & deg$P.Value < 0.05, ] # |log2FC| ≥ 1, FDR < 5%

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
        file = paste0("binary_deg_", subtype, "_", bulk_i, ".csv"),
        row.names = TRUE
      )
    }
  )

  gc()

  # ? Find markers for survival phenotype
  future::plan(future::multicore, workers = 8)
  options(future.globals.maxSize = 1024^3 * 3)
  cli::cli_h2("Finding marker genes for survival phenotype")

  purrr::walk(
    s_bulks,
    function(bulk_i) {
      bulk <- qs::qread(
        file.path(data_path, paste0("brca_bulkdata_", bulk_i, ".qs")),
        nthreads = 4L
      )

      bulk_genes <- rownames(bulk)
      cm_genes <- intersect(sc_genes, bulk_genes)

      cli::cli_alert_info(
        "{bulk_i} has {.val {length(cm_genes)}} common genes with single cell data"
      )

      bulk <- bulk[cm_genes, ]

      if (bulk_i == "GSE42568") {
        surv <- qs::qread(
          file.path(data_path, paste0("brca_pheno_", bulk_i, ".qs")),
          nthreads = 4
        ) %>%
          select(
            "overall survival time_days:ch1",
            "overall survival event:ch1"
          ) %>%
          filter(`overall survival event:ch1` != "NA") %>% # cannot be changed to `!is.na()`
          rename("time" := 1, "status" := 2) %>%
          mutate_all(~ as.numeric(.))
      } else if (bulk_i == "GSE162228") {
        surv <- qs::qread(
          file.path(data_path, paste0("brca_pheno_", bulk_i, ".qs")),
          nthreads = 4
        ) %>%
          select("overall survival (years):ch1", "characteristics_ch1.5") %>%
          rename("time" := 1) %>%
          mutate(
            status = case_when(
              characteristics_ch1.5 == "alive: Death" ~ 1,
              characteristics_ch1.5 == "alive: Alive" ~ 0
            )
          ) %>%
          select(-"characteristics_ch1.5") %>%
          mutate_all(~ as.numeric(.))
      } else if (bulk_i == "TCGA") {
        surv <- qs::qread(file.path(data_path, "brca_surv_TCGA.qs"))
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
      dt[is.na(gene), gene := paste0("unknown_", .I)] # remove NA
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
            HR_low = exp(
              stats::coef(fit) - 1.96 * sqrt(stats::vcov(fit)[1, 1])
            ),
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
        dplyr::filter(sig)

      if (nrow(filtered_cox_results) == 0) {
        cli::cli_h2(cli::col_red(
          "{bulk_i} has no significant genes in survival"
        ))
      } else if (nrow(filtered_cox_results) < 50) {
        cli::cli_h2(cli::col_red(
          "{bulk_i} has less than 20 significant genes in survival"
        ))
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
        file = paste0("survival_deg_", subtype, "_", bulk_i, ".csv")
      )
    }
  )
})
