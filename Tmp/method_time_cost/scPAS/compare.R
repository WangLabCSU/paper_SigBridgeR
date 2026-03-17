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
  file.path(data_root, "ov_bulkdata_GSE140082.qs"),
  nthreads = 2L
)
# GSM4153781 has NA values, find median and replace NA
median = median(bulk[, "GSM4153781"], na.rm = TRUE)
na_indices <- which(is.na(bulk[, "GSM4153781"]))
bulk[na_indices, "GSM4153781"] <- median

pheno <- qs::qread(
  file.path(data_root, "ov_pheno_GSE140082.qs"),
  nthreads = 2L
)

surv <- pheno %>%
  select("final_ostm:ch1", "final_osid:ch1") %>%
  rename("time" = 1, "status" = 2) %>%
  mutate_all(~ as.numeric(.)) %>%
  filter(rownames(.) %in% colnames(bulk))

intersect(rownames(seurat), rownames(bulk)) |> length()
# 11620

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

# seurat_100k <- qs::qread(
#   "/home/data/sigbridger/method_time_cost/seurat_100k.qs",
#   nthreads = 4L
# )
# seurat_100k <- SeuratObject::JoinLayers(seurat_100k)

#  ----------------- Benchmark -----------------------

survival_tweaked <- press(
  n_cells = c(
    1e3,
    5e3,
    1e4,
    5e4
    # , 1e5
  ),
  {
    seurat_obj <- switch(
      as.character(n_cells),
      "1000" = seurat_1k,
      "5000" = seurat_5k,
      "10000" = seurat_10k,
      "50000" = seurat_50k
      # ,      "100000" = seurat_100k
    )

    mark(
      Screen(
        matched_bulk = bulk,
        sc_data = seurat_obj,
        phenotype = surv,
        phenotype_class = "survival",
        screen_method = "scPAS"
      ),
      check = FALSE,
      iterations = 3
    )
  }
) %>%
  mutate(n_cells = as.factor(n_cells))

qs::qsave(
  survival_tweaked,
  file.path(output_dir, "scPAS_survival_bench.qs"),
  nthreads = 8L
)

gc()

cli::cli_h1("scPAS survival benchmark done")
