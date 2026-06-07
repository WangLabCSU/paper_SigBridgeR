# ! TCGA_LUAD
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

distance_choices <- c(
  "cosine",
  "pearson",
  "spearman",
  "kendall",
  "euclidean",
  "maximum"
)

# * random search, 50 times
set.seed(123)
arg_samples <- data.frame(
  distance = sample(distance_choices, 50, replace = TRUE), # 第1维
  nPerm = sample(seq(200, 2500, 100), 50, replace = TRUE),
  log2FC = sample(seq(1.2, 3, 0.01), 50, replace = TRUE)
) %>%
  dplyr::add_row(distance = "cosine", nPerm = 1000L, log2FC = 1L) # default parameters


options(future.globals.maxSize = 30 * 1024^3)
future::plan(future.mirai::mirai_multisession(workers = 4L))
SigBridgeR::setThreads(4L)

# ! To avoid recomputing, file cache is used
if (!dir.exists("stats/pipet1")) {
  dir.create("stats/pipet1", recursive = TRUE)
}


res_list <- lapply(
  seq_len(nrow(arg_samples)),
  function(i) {
    cli::cli_h1("{i} / {nrow(arg_samples)}")

    # ! load cache if exists
    cache_save_path <- file.path("stats/pipet1", glue::glue("process_{i}.csv"))
    if (file.exists(cache_save_path)) {
      cli::cli_alert("cache found, loading...")
      cache <- data.table::fread(cache_save_path)
      return(cache)
    }

    result <- suppressWarnings(SigBridgeR::Screen(
      matched_bulk = bulk,
      sc_data = sc_data,
      phenotype = pheno_bi,
      label_type = glue::glue("process_{i}"),
      phenotype_class = "binary",
      screen_method = "PIPET",
      distance = arg_samples$distance[i], # select_alpha will be used
      nPerm = as.integer(arg_samples$nPerm[i]),
      log2FC = arg_samples$log2FC[i],
      verbose = FALSE,
      parallel = TRUE
    ))

    data <- data.frame(
      pos_cell = (result$scRNA_data$PIPET == "Positive")
    )
    colnames(data) <- glue::glue("process_{i}")

    # ! save cache
    data.table::fwrite(data, cache_save_path)

    # 返回包含索引和结果的数据框
    return(data)
  }
)


# 合并所有结果
all_results <- do.call(cbind, res_list)
rownames(all_results) = colnames(sc_data) # each cell is a row

data.table::fwrite(
  all_results,
  file = "stats/pipet_label_mat1.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("(1) pipet random search completed."))

# ! TCGA_LUAD
