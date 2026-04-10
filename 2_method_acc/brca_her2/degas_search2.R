# ! GSE42568

setwd(file.path(usethis::proj_path(), "2_method_acc/brca_her2"))
library(dplyr)
library(SigBridgeR)


# * Load Data
data_dir <- "/home/data/sigbridger/benchmark_data/brca"

sc_data <- qs::qread(file.path(data_dir, "seurat_her2.qs"), nthreads = 8L)

bulk <- qs::qread(
  file.path(data_dir, "brca_bulkdata_GSE42568.qs"),
  nthreads = 2L
)

cli::cli_alert_info("bulk data loaded: dim = ({.val {dim(bulk)}})")

pheno <- qs::qread(file.path(data_dir, "brca_pheno_GSE42568.qs"))

cm_samples <- intersect(rownames(pheno), colnames(bulk))


pheno_bi <- setNames(
  ifelse(pheno$`tissue:ch1` == "breast cancer", 1L, 0L),
  cm_samples
)
bulk <- bulk[, names(pheno_bi)]


cli::cli_alert_info("pheno data loaded: 1~tumor, 0~normal")
table(pheno_bi)

if (!all(names(pheno_bi) == colnames(bulk))) {
  stop("pheno_bi and bulk not match")
}

seurat_tumor <- readRDS(
  "/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_HER2Tum.rds"
)
tumor_cells <- rownames(seurat_tumor@meta.data)
benchmark_label <- colnames(sc_data) %in% tumor_cells

# * Screen

SigBridgeR::setThreads(
  8L,
  tf_config = list(
    xla = TRUE,
    intra_op = 8L,
    inter_op = 8L
  )
)


# * random search, 50 times
set.seed(123)
arg_samples <- data.frame(
  arch = sample(c("DenseNet", "Standard"), 50, replace = TRUE),
  ff_depth = sample(2:10, 50, replace = TRUE),
  bag_depth = sample(3:10, 50, replace = TRUE)
) %>%
  dplyr::add_row(arch = "DenseNet", ff_depth = 3, bag_depth = 5)

# * run DEGAS
res_list <- lapply(
  seq_len(nrow(arg_samples)),
  function(i) {
    cli::cli_h1("{i} / {nrow(arg_samples)}")

    tryCatch(
      {
        arch <- arg_samples[i, "arch"][[1]]
        ff_depth <- arg_samples[i, "ff_depth"][[1]]
        bag_depth <- arg_samples[i, "bag_depth"][[1]]

        result <- Screen(
          matched_bulk = bulk,
          sc_data = sc_data,
          phenotype = pheno_bi,
          label_type = glue::glue("process_{i}"),
          phenotype_class = "binary",
          screen_method = "DEGAS",
          degas_params = list(
            DEGAS.architecture = arch,
            DEGAS.ff_depth = ff_depth,
            DEGAS.bag_depth = bag_depth
          )
        )

        data <- data.frame(
          pos_cell = (result$scRNA_data$DEGAS == "Positive")
        )
        colnames(data) <- glue::glue("process_{i}")
        gc()

        # 返回包含索引和结果的数据框
        return(data)
      },
      error = function(e) {
        cli::cli_alert_danger("ERROR {i}: {e$message}")
        data <- data.frame(
          pos_cell = FALSE
        )
        colnames(data) = glue::glue("process_{i}")
        return(data)
      }
    )
  }
)

gc()
all_results <- do.call(cbind, res_list)
rownames(all_results) <- colnames(sc_data)

data.table::fwrite(
  all_results,
  file = "stats/degas_label_mat2_part1.csv",
  row.names = TRUE
)

# -------------------------------------------------------------------------------------------------

arg_samples2 <- data.frame(
  lamb1 = sample(2:10, 50, replace = TRUE),
  lamb2 = sample(2:10, 50, replace = TRUE),
  lamb3 = sample(2:10, 50, replace = TRUE)
) %>%
  dplyr::add_row(lamb1 = 3, lamb2 = 3, lamb3 = 3)

# * run DEGAS
res_list <- lapply(
  seq_len(nrow(arg_samples2)),
  function(i) {
    cli::cli_h1("{i} / {nrow(arg_samples2)}")

    tryCatch(
      {
        lamb1 <- arg_samples2[i, "lamb1"][[1]]
        lamb2 <- arg_samples2[i, "lamb2"][[1]]
        lamb3 <- arg_samples2[i, "lamb3"][[1]]

        result <- Screen(
          matched_bulk = bulk,
          sc_data = sc_data,
          phenotype = pheno_bi,
          label_type = glue::glue("process_{i}"),
          phenotype_class = "binary",
          screen_method = "DEGAS",
          degas_params = list(
            DEGAS.lambda1 = lamb1,
            DEGAS.lambda2 = lamb2,
            DEGAS.lambda3 = lamb3
          )
        )

        pos <- (result$scRNA_data$DEGAS == "Positive")

        data <- data.frame(
          pos_cell = pos
        )
        colnames(data) = glue::glue("process_{i}")
        gc()

        # 返回包含索引和结果的数据框
        return(data)
      },
      error = function(e) {
        cli::cli_alert_danger("ERROR {i}: {e$message}")
        data <- data.frame(
          pos_cell = FALSE
        )
        colnames(data) = glue::glue("process_{i}")
        return(data)
      }
    )
  }
)

gc()
all_results <- do.call(cbind, res_list)
rownames(all_results) <- colnames(sc_data)

data.table::fwrite(
  all_results,
  file = "stats/degas_label_mat2_part2.csv",
  row.names = TRUE
)

# -------------------------------------------------------------------------------------------------

arg_samples3 <- data.frame(
  scbatch_sz = sample(50:500, 50, replace = TRUE),
  patbatch_sz = sample(25:100, 50, replace = TRUE),
  hidden_feats = sample(25:100, 50, replace = TRUE),
  do_prc = sample(seq(0.1, 0.9, 0.1), 50, replace = TRUE)
) %>%
  dplyr::add_row(
    scbatch_sz = 200,
    patbatch_sz = 50,
    hidden_feats = 50,
    do_prc = 0.5
  )

res_list <- lapply(
  seq_len(nrow(arg_samples3)),
  function(i) {
    cli::cli_h1("{i} / {nrow(arg_samples3)}")
    tryCatch(
      {
        scbatch_sz <- arg_samples3[i, "scbatch_sz"][[1]]
        patbatch_sz <- arg_samples3[i, "patbatch_sz"][[1]]
        hidden_feats <- arg_samples3[i, "hidden_feats"][[1]]
        do_prc <- arg_samples3[i, "do_prc"][[1]]

        result <- Screen(
          matched_bulk = bulk,
          sc_data = sc_data,
          phenotype = pheno_bi,
          label_type = glue::glue("process_{i}"),
          phenotype_class = "binary",
          screen_method = "DEGAS",
          degas_params = list(
            DEGAS.scbatch_sz = scbatch_sz,
            DEGAS.patbatch_sz = patbatch_sz,
            DEGAS.hidden_feats = hidden_feats,
            DEGAS.do_prc = do_prc
          )
        )

        pos <- (result$scRNA_data$DEGAS == "Positive")

        data <- data.frame(
          pos_cell = pos
        )
        colnames(data) = glue::glue("process_{i}")
        gc()

        # 返回包含索引和结果的数据框
        return(data)
      },
      error = function(e) {
        cli::cli_alert_danger("ERROR {i}: {e$message}")
        data <- data.frame(
          pos_cell = FALSE
        )
        colnames(data) = glue::glue("process_{i}")
        return(data)
      }
    )
  }
)

gc()
all_results <- do.call(cbind, res_list)
rownames(all_results) = colnames(sc_data)

data.table::fwrite(
  all_results,
  file = "stats/degas_label_mat2_part3.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("(2)DEGAS random search completed."))
