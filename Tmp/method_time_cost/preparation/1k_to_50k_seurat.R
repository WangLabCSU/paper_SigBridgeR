data_root <- "/home/data/sigbridger/benchmark_data/ov"

seurat <- qs::qread(
  file.path(data_root, "hgsoc_GSE165897_seurat.qs"),
  nthreads = 8L
)

set.seed(123L)

seurat_1k <- subset(seurat, cells = sample(colnames(seurat), size = 1e3))

seurat_5k <- subset(seurat, cells = sample(colnames(seurat), size = 5e3))

seurat_10k <- subset(seurat, cells = sample(colnames(seurat), size = 1e4))

seurat_50k <- subset(seurat, cells = sample(colnames(seurat), size = 5e4))

out_dir <- "/home/data/sigbridger/method_time_cost"


purrr::walk2(
  list(seurat_1k, seurat_5k, seurat_10k, seurat_50k),
  list("seurat_1k", "seurat_5k", "seurat_10k", "seurat_50k"),
  ~ qs::qsave(
    .x,
    file.path(out_dir, paste0(.y, ".qs")),
    nthreads = 8L
  )
)
