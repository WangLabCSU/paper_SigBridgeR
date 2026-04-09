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
set.seed(12345)

arg_sample <- data.frame(
  n_pc = sample(seq(10L, 100L, 10L), 50, replace = TRUE),
  resolution = sample(seq(0.5, 5L, 0.5), 50, replace = TRUE),
  ela_net_alpha = sample(seq(0.1, 1L, 0.1), 50, replace = TRUE)
) %>%
  dplyr::add_row(n_pc = 60L, resolution = 2L, ela_net_alpha = 0.4) # default

res_list <- lapply(
  seq_len(nrow(arg_sample)),
  function(i) {
    n_pc_i <- arg_samples$n_pc[i]
    resolution_i <- arg_samples$resolution[i]
    ela_net_alpha_i <- arg_samples$ela_net_alpha[i]

    tryCatch(
      {
        result <- SigBridgeR::Screen(
          bulk,
          sc_data,
          pheno_bi,
          screen_method = "SCIPAC",
          label_type = glue::glue("process_{i}"),
          phenotype_class = "binary",
          n_pc = n_pc_i,
          resolution = resolution_i,
          ela_net_alpha = ela_net_alpha_i
        )
        pos_cell <- (result$scRNA_data$SCIPAC == "Positive")

        data <- data.frame(
          pos_cell = pos_cell
        )

        colnames(data) <- glue::glue("process_{i}")
        gc(verbose = FALSE)

        # 返回包含索引和结果的数据框
        return(data)
      },
      error = function(e) {
        cli::cli_alert_warning(c(
          "x" = "SCIPAC result is not complete, maybe not suitable for this parameter pair, using all FALSE"
        ))
        data <- data.frame(pos_cell = rep(FALSE, ncol(sc_data)))
        colnames(data) <- glue::glue("process_{i}")
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
  file = "stats/scipac_label_mat1.csv",
  row.names = TRUE
)

cli::cli_alert_success(crayon::green("(1) scipac random search completed."))

# ! TCGA_BRCA
