# ! TCGA_LUAD

library(SigBridgeR)
library(dplyr)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  file.path(usethis::proj_path(), "2_method_acc/lung")
)
data_path <- "/home/data/sigbridger/benchmark_data/lung"

# * load data
sc_data <- qs::qread(file.path(data_path, "luad_GSE123902_seurat.qs"))

bulk <- qs::qread(
  file.path(data_path, "TCGA_LUAD_bulkdata.qs")
)
bulk <- log2(bulk + 1)

pheno <- qs::qread(file.path(data_path, "TCGA_LUAD_pheno.qs"))

pheno_bi <- mutate(pheno, sample_type = substr(pheno$sample, 14, 15)) %>%
  select(sample, sample_type) %>%
  filter(sample_type %in% c("01", "11")) %>%
  mutate(sample_type = as.integer(sample_type == "01"))
pheno_bi <- setNames(pheno_bi$sample_type, pheno_bi$sample)

bulk <- bulk[, names(pheno_bi)]

if (!all(colnames(bulk) == names(pheno_bi))) {
  stop("bulk and pheno_bi not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}
if (anyNA(pheno_bi)) {
  stop("pheno_bi has NA")
}

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
if (!file.exists("stats/degas_label_mat1_part1.csv")) {
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
            label_type = "DEGAS_arg_one_",
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
    file = "stats/degas_label_mat1_part1.csv",
    row.names = TRUE
  )
} else {
  cli::cli_alert_info("Found degas_label_mat1_part1.csv, skip")
}


# -------------------------------------------------------------------------------------------------

arg_samples2 <- data.frame(
  lamb1 = sample(2:10, 50, replace = TRUE),
  lamb2 = sample(2:10, 50, replace = TRUE),
  lamb3 = sample(2:10, 50, replace = TRUE)
) %>%
  dplyr::add_row(lamb1 = 3, lamb2 = 3, lamb3 = 3)

# * run DEGAS
if (!file.exists("stats/degas_label_mat1_part2.csv")) {
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
            label_type = "DEGAS_arg_two_",
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
    file = "stats/degas_label_mat1_part2.csv",
    row.names = TRUE
  )
} else {
  cli::cli_alert_info("Found degas_label_mat1_part2.csv, skip")
}


# -------------------------------------------------------------------------------------------------

arg_samples3 <- data.frame(
  scbatch_sz = sample(seq(50, 500, 10), 50, replace = TRUE),
  patbatch_sz = sample(seq(25, 100, 5), 50, replace = TRUE),
  hidden_feats = sample(seq(25, 100, 5), 50, replace = TRUE),
  do_prc = sample(seq(0.1, 0.9, 0.1), 50, replace = TRUE)
) %>%
  dplyr::add_row(
    scbatch_sz = 200,
    patbatch_sz = 50,
    hidden_feats = 50,
    do_prc = 0.5
  )

# ! To avoid recomputing, file cache is used
if (!dir.exists("stats/degas/part3")) {
  dir.create("stats/degas/part3", recursive = TRUE)
}


if (!file.exists("stats/degas_label_mat1_part3.csv")) {
  res_list <- lapply(
    seq_len(nrow(arg_samples3)),
    function(i) {
      cli::cli_h1("{i} / {nrow(arg_samples3)}")

      # ! load cache if exists
      cache_save_path <- file.path(
        "stats/degas/part3",
        glue::glue("process_{i}.csv")
      )
      if (file.exists(cache_save_path)) {
        cli::cli_alert("cache found, loading...")
        cache <- data.table::fread(cache_save_path)
        return(cache)
      }

      scbatch_sz <- arg_samples3[i, "scbatch_sz"][[1]]
      patbatch_sz <- arg_samples3[i, "patbatch_sz"][[1]]
      hidden_feats <- arg_samples3[i, "hidden_feats"][[1]]
      do_prc <- round(arg_samples3[i, "do_prc"][[1]], 3)

      result <- Screen(
        matched_bulk = bulk,
        sc_data = sc_data,
        phenotype = pheno_bi,
        label_type = "DEGAS_arg_three_",
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

      # ! save cache
      data.table::fwrite(data, cache_save_path)

      # 返回包含索引和结果的数据框
      return(data)
    }
  )
  gc()
  all_results <- do.call(cbind, res_list)
  rownames(all_results) = colnames(sc_data)

  data.table::fwrite(
    all_results,
    file = "stats/degas_label_mat1_part3.csv",
    row.names = TRUE
  )
} else {
  cli::cli_alert_info("Found degas_label_mat1_part3.csv, skip")
}

# ! TCGA_LUAD
