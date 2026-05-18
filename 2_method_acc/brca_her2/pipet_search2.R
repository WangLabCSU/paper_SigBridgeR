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
  as.integer(pheno$`tissue:ch1` == "breast cancer"),
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
# SigBridgeR::setThreads(4L)

# ! To avoid recomputing, file cache is used
if (!dir.exists("stats/pipet2")) {
  dir.create("stats/pipet2", recursive = TRUE)
}

res_list <- lapply(
  seq_len(nrow(arg_samples)),
  function(i) {
    cli::cli_h1("{i} / {nrow(arg_samples)}")

    # ! load cache if exists
    cache_save_path <- file.path("stats/pipet2", glue::glue("process_{i}.csv"))
    if (file.exists(cache_save_path)) {
      cli::cli_alert("cache found, loading...")
      cache <- data.table::fread(cache_save_path)
      return(cache)
    }

    result <- suppressWarnings(Screen(
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
  file = "stats/pipet_label_mat2.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("(2) pipet random search completed."))

# ! GSE42568
