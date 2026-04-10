# ! TCGA_BRCA

setwd(file.path(usethis::proj_path(), "2_method_acc/brca_her2"))
library(dplyr)

# * Load Data
data_dir <- "/home/data/sigbridger/benchmark_data/brca"

sc_data <- qs::qread(file.path(data_dir, "seurat_her2.qs"), nthreads = 8L)

bulk <- qs::qread(file.path(data_dir, "brca_bulkdata_TCGA.qs"), nthreads = 2L)
bulk <- log2(bulk + 1)
cli::cli_alert_info("bulk data loaded: dim = ({.val {dim(bulk)}})")

pheno <- qs::qread(file.path(data_dir, "brca_pheno_TCGA.qs"))

cm_samples <- intersect(pheno$sample, colnames(bulk))

pheno_bi <- pheno %>%
  mutate(sample_type = substr(sample, 14, 15)) %>%
  filter(sample_type %in% c("01", "11"), sample %in% cm_samples) %>%
  mutate(sample_type = as.integer(sample_type == "01")) %>%
  {
    setNames(.$sample_type, .$sample)
  }

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
# ! scPP random search, 50 times

# * random search, 51 times
set.seed(12345)

options(future.globals.maxSize = 20 * 1024^3)
future::plan(future::multicore, workers = 4L)

# ! To avoid recomputing, file cache is used
if (!dir.exists("stats/scpp1")) {
  dir.create("stats/scpp1", recursive = TRUE)
}


arg_samples <- data.frame(
  prob = sample(seq(0, 0.5, 0.01), 50, replace = TRUE),
  Log2FC_cutoff = round(runif(50), 3)
) %>%
  dplyr::add_row(prob = 0.2, Log2FC_cutoff = 0.585) # default

res_list <- lapply(
  seq_len(nrow(arg_samples)),
  function(i) {
    cli::cli_h1("{i} / {nrow(arg_samples)}")

    # ! load cache if exists
    cache_save_path <- file.path("stats/scpp1", glue::glue("process_{i}.csv"))
    if (file.exists(cache_save_path)) {
      cli::cli_alert("cache found, loading...")
      cache <- data.table::fread(cache_save_path)
      return(cache)
    }

    prob_i <- arg_samples$prob[i]
    Log2FC_cutoff_i <- arg_samples$Log2FC_cutoff[i]

    tryCatch(
      {
        scpp_result <- SigBridgeR::Screen(
          bulk,
          sc_data,
          pheno_bi,
          screen_method = "scPP",
          label_type = glue::glue("process_{i}"),
          phenotype_class = "binary",
          ref_group = 0L,
          Log2FC_cutoff = Log2FC_cutoff_i,
          probs = prob_i,
          parallel = TRUE,
          assay = "RNA"
        )
        pos_cell <- (scpp_result$scRNA_data$scPP == "Positive")

        data <- data.frame(
          pos_cell = pos_cell
        )

        colnames(data) <- glue::glue("process_{i}")
        gc(verbose = FALSE)

        # ! save cache
        data.table::fwrite(data, cache_save_path)

        # 返回包含索引和结果的数据框
        return(data)
      },
      error = function(e) {
        cli::cli_alert_warning(c(
          "x" = "scPP result is not complete, maybe not suitable for this parameter pair, using all FALSE.",
          ">" = "prob = {prob_i}, Log2FC_cutoff = {Log2FC_cutoff_i}, process_{i}"
        ))
        data <- data.frame(pos_cell = rep(FALSE, ncol(sc_data)))
        colnames(data) <- glue::glue("process_{i}")

        # ! save cache
        data.table::fwrite(data, cache_save_path)

        return(data)
      }
    )
  }
)


# *visualize
gc()
all_results <- do.call(cbind, res_list)
rownames(all_results) = colnames(sc_data)

data.table::fwrite(
  all_results,
  file = "stats/scpp_label_mat1.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("(1) scpp random search completed."))

# ! TCGA_BRCA
