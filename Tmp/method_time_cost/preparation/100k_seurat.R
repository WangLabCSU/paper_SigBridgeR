data_path <- "/home/data/data-resource/single-cell/NC-Zenodo10651059/atlas_dataset"
# * get all ov cancer
all_sc_paths <- list.files(
  data_path,
  pattern = "ov_hgs.*\\.h5ad$",
  all.files = TRUE,
  full.names = TRUE
)
all_sc_paths
# [1] "/home/data/data-resource/single-cell/NC-Zenodo10651059/atlas_dataset/ov_hgs_GSE154600.h5ad"
# [2] "/home/data/data-resource/single-cell/NC-Zenodo10651059/atlas_dataset/ov_hgs_GSE158937.h5ad"
# [3] "/home/data/data-resource/single-cell/NC-Zenodo10651059/atlas_dataset/ov_hgstoc_pan_blueprint_ovary.h5ad"

all_sc_paths <- all_sc_paths[1:2]

all_sc_names <- basename(all_sc_paths) |> tools::file_path_sans_ext()

purrr::walk2(all_sc_paths, all_sc_names, \(path, name) {
  adata <- anndataR::read_h5ad(path)
  X <- Matrix::t(adata$X)
  counts <- exp(X) - 1

  seurat <- SeuratObject::CreateSeuratObject(
    counts,
    meta = adata$obs,
    min.cells = 400
  )
  assay_data <- SeuratObject::CreateSeuratObject(counts = X, min.cells = 400) # data
  seurat@assays$RNA$data <- assay_data@assays$RNA$counts # data

  seurat <- Seurat::FindVariableFeatures(seurat)
  var_features <- SeuratObject::VariableFeatures(seurat)
  seurat <- Seurat::ScaleData(seurat)
  seurat <- Seurat::RunPCA(seurat, features = var_features)
  #   seurat <- Seurat::FindNeighbors(seurat, dims = 1:10)
  assign(name, seurat, envir = .GlobalEnv)
})

# -----------------------------------------------------------------------------------------------

out_dir <- "/home/data/sigbridger/method_time_cost"

# * We will use the first GSE154600, GSE158937 and GSE165897

# qs::qsave(ov_hgs_GSE165897, file.path(out_dir, "gse165897.qs"), nthreads = 4L)
qs::qsave(ov_hgs_GSE154600, file.path(out_dir, "gse154600.qs"), nthreads = 4L)
qs::qsave(ov_hgs_GSE158937, file.path(out_dir, "gse158937.qs"), nthreads = 4L)

gse165897 <- qs::qread(
  "/home/data/sigbridger/benchmark_data/ov/hgsoc_GSE165897_seurat.qs",
  nthreads = 4L
)
gse154600 <- qs::qread(file.path(out_dir, "gse154600.qs"), nthreads = 4L)
gse158937 <- qs::qread(file.path(out_dir, "gse158937.qs"), nthreads = 4L)

seurat_100k <- SCIntegrate.Seurat(
  gse165897,
  gse154600,
  gse158937,
  method = Seurat::HarmonyIntegration,
  new.reduction = "harmony",
  dims = 1:10, # passed to FindNeighbors, RunUMAP, etc.
  resolution = 0.5
)

# * OK
get_names_4_ids(gse165897, gse154600, gse158937)

merged <- merge(
  x = gse165897,
  y = c(gse154600, gse158937),
  add.cell.ids = get_names_4_ids(gse165897, gse154600, gse158937),
  project = "Integrated Seurat"
)


merged <- Seurat::FindVariableFeatures(merged)
merged <- Seurat::ScaleData(merged)
# variable_features <- Seurat::VariableFeatures(merged)
merged <- Seurat::RunPCA(merged)

merged <- rlang::exec(
  Seurat::IntegrateLayers,
  object = merged,
  method = Seurat::HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = TRUE
)
merged <- Seurat::FindNeighbors(
  merged,
  reduction = "harmony",
  dims = 1:10
)
merged <- Seurat::FindClusters(merged, resolution = 0.6)

qs::qsave(merged, file.path(out_dir, "seurat_100k.qs"), nthreads = 4L)
merged <- qs::qread(file.path(out_dir, "seurat_100k.qs"), nthreads = 4L)

merged <- SeuratObject::JoinLayers(merged)

qs::qsave(merged, file.path(out_dir, "seurat_100k_joined.qs"), nthreads = 4L)
