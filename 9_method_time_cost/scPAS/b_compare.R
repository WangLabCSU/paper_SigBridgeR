library(bench)
library(dplyr)
library(SigBridgeR)

#  ----------------- Load Data -----------------------

data_root <- "/home/data/sigbridger/benchmark_data/ov"

seurat <- qs::qread(
  file.path(data_root, "hgsoc_GSE165897_seurat.qs"),
  nthreads = 8L
)

dim(seurat)
# [1] 15201 50865

bulk <- qs::qread(
  file.path(data_root, "ov_bulkdata_GSE9891.qs"),
  nthreads = 2L
)

pheno <- qs::qread(
  file.path(data_root, "ov_pheno_GSE9891.qs"),
  nthreads = 2L
)

pheno_bi <- setNames(
  case_when(
    pheno$characteristics_ch1.1 == "Type : LMP" ~ 0,
    pheno$characteristics_ch1.1 == "Type : Malignant" ~ 1
  ),
  pheno$geo_accession
)

bulk_bi <- bulk[, names(pheno_bi)]

# -------------------------------- Output ------------------------------------------------

output_dir <- "/home/data/sigbridger/method_time_cost/scPAS"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#  ----------------- Subset Single Cell Data -----------------------

set.seed(123L)

seurat_1k <- subset(seurat, cells = sample(colnames(seurat), size = 1e3))

seurat_5k <- subset(seurat, cells = sample(colnames(seurat), size = 5e3))

seurat_10k <- subset(seurat, cells = sample(colnames(seurat), size = 1e4))

seurat_50k <- subset(seurat, cells = sample(colnames(seurat), size = 5e4))

seurat_100k <- qs::qread(
  "/home/data/sigbridger/method_time_cost/seurat_100k.qs",
  nthreads = 4L
)
seurat_100k <- SeuratObject::JoinLayers(seurat_100k)

#  ----------------- Benchmark -----------------------

binary_tweaked <- press(
  n_cells = c(1e3, 5e3, 1e4, 5e4, 1e5),
  {
    seurat_obj <- switch(
      as.character(n_cells),
      "1000" = seurat_1k,
      "5000" = seurat_5k,
      "10000" = seurat_10k,
      "50000" = seurat_50k,
      "100000" = seurat_100k
    )

    mark(
      Screen(
        matched_bulk = bulk_bi,
        sc_data = seurat_obj,
        phenotype = pheno_bi,
        phenotype_class = "binary",
        screen_method = "scPAS"
      ),
      check = FALSE,
      iterations = 3
    )
  }
) %>%
  mutate(n_cells = as.factor(n_cells))

gc()

qs::qsave(
  binary_tweaked,
  file.path(output_dir, "scPAS_binary_bench.qs"),
  nthreads = 8L
)


cli::cli_h1("Finished benchmarking scPAS!")
