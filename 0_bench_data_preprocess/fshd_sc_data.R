# Read raw scRNA-seq counts of FSHD dataset GSE122873
# Each GSM file is a gene (Ensembl ID, rows) x cell (barcode, cols) raw count matrix
# Samples: FSHD1.1, FSHD1.2, FSHD2.1, FSHD2.2, CTRL.1, CTRL.2

fshd_dir <- "/home/data/sigbridger/benchmark_data/fshd/GSE122873/GSE122873_RAW"
fshd_files <- list.files(fshd_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)

read_raw_counts <- function(file) {
  dt <- data.table::fread(file)
  genes <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- genes
  Matrix::Matrix(mat, sparse = TRUE)
}

fshd_counts_list <- lapply(fshd_files, function(file) {
  # sample name, e.g. "FSHD1.1" from "GSM3487556_FSHD1.1.txt.gz"
  sample <- sub("^GSM\\d+_(.+)\\.txt\\.gz$", "\\1", basename(file))
  counts <- read_raw_counts(file)
  # barcode collision across samples: append sample as suffix
  colnames(counts) <- paste0(colnames(counts), "_", sample)
  counts
})
names(fshd_counts_list) <- sub(
  "^GSM\\d+_(.+)\\.txt\\.gz$",
  "\\1",
  basename(fshd_files)
)

# metadata: sample / condition per cell
fshd_meta <- data.frame(
  barcode = unlist(lapply(names(fshd_counts_list), function(s) {
    colnames(fshd_counts_list[[s]])
  })),
  sample = purrr::imap(fshd_counts_list, \(x, name) {
    rep(name, ncol(x))
  }) %>%
    unlist(),
  row.names = NULL
)
fshd_meta$condition <- sub("\\.\\d+$", "", fshd_meta$sample) # FSHD1 / FSHD2 / CTRL
rownames(fshd_meta) <- fshd_meta$barcode

# merged raw counts (all genes shared across samples)
common_genes <- Reduce(intersect, lapply(fshd_counts_list, rownames))
stringr::str_detect(common_genes, "\\.") %>% any()
fshd_counts <- do.call(
  cbind,
  lapply(fshd_counts_list, function(x) x[common_genes, ])
)

# sample-level metadata from GEO (disease state, sex, ...)
fshd_meta_plus <- Biobase::pData(geokit::geo_matrix(
  "GSE122873",
  odir = "/home/data/sigbridger/benchmark_data/fshd/GSE122873",
  add_gpl = TRUE
))

# merge sample-level info into cell-level meta (match by sample name, e.g. "FSHD1.1")
fshd_meta <- merge(
  fshd_meta,
  fshd_meta_plus,
  by.x = "sample",
  by.y = "title",
  all.x = TRUE,
  sort = FALSE
)
rownames(fshd_meta) <- fshd_meta$barcode
fshd_meta <- fshd_meta[colnames(fshd_counts), ] # keep cell order consistent

# Seurat object with raw counts
fshd_seurat <- Seurat::CreateSeuratObject(
  counts = fshd_counts,
  meta.data = fshd_meta,
  project = "FSHD_GSE122873"
)

# ---- Standard Seurat (v5) preprocessing -----------------------------------

# map Ensembl IDs to gene symbols (needed to detect mitochondrial genes)
ens2sym <- IDConverter::convert_hm_genes(rownames(fshd_seurat), "ensembl")
cli::cli_alert_warning("Found {sum(is.na(ens2sym))} NA value{?s}")

mt_genes_found <- rownames(fshd_seurat)[
  !is.na(ens2sym) & grepl("^MT-", ens2sym)
]
cli::cli_alert_info("Found {length(mt_genes_found)} mitochondrial genes")

# QC metrics overview BEFORE filtering
qc_plot_before <- Seurat::VlnPlot(
  fshd_seurat,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "sample",
  ncol = 2,
  pt.size = 0
)

# QC filtering with MAD-based adaptive thresholds (OSCA convention), per sample:
# outlier = beyond median +/- 3 * MAD (log-scale for library size / nFeature)
qc <- fshd_seurat[[]]
qc$log10_nFeature <- log10(qc$nFeature_RNA)
qc$log10_nCount <- log10(qc$nCount_RNA)

low_lib <- scuttle::isOutlier(
  qc$log10_nCount,
  nmads = 3,
  type = "lower",
  batch = qc$sample
)
low_feat <- scuttle::isOutlier(
  qc$log10_nFeature,
  nmads = 3,
  type = "lower",
  batch = qc$sample
)
qc$qc_pass <- !(low_lib | low_feat)

# per-sample adaptive threshold summary + removal stats
qc_summary <- lapply(split(qc, qc$sample), function(d) {
  med_nf <- median(d$log10_nFeature)
  mad_nf <- mad(d$log10_nFeature)
  med_mt <- median(d$percent.mt)
  mad_mt <- mad(d$percent.mt)
  data.frame(
    sample = d$sample[1],
    n_cells = nrow(d),
    nFeature_lower = round(10^(med_nf - 3 * mad_nf)),
    n_removed = sum(!d$qc_pass),
    pct_removed = round(mean(!d$qc_pass) * 100, 2)
  )
}) %>%
  dplyr::bind_rows()

print(qc_summary)
data.table::fwrite(
  qc_summary,
  "/home/data/sigbridger/benchmark_data/fshd/GSE122873/fshd_sc_GSE122873_qc_threshold_summary.csv"
)

# visualize adaptive thresholds on top of QC metrics
qc_scatter <- ggplot2::ggplot(
  qc,
  ggplot2::aes(nCount_RNA, nFeature_RNA, color = qc_pass)
) +
  ggplot2::geom_point(size = 0.4, alpha = 0.5) +
  ggplot2::facet_wrap(~sample, scales = "free") +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10() +
  ggplot2::theme_bw()


fshd_seurat <- subset(fshd_seurat, cells = rownames(qc)[qc$qc_pass])
fshd_seurat <- subset(
  fshd_seurat,
  features = rownames(fshd_seurat)[
    Matrix::rowSums(Seurat::GetAssayData(fshd_seurat, layer = "counts") > 0) >=
      3
  ]
)


# find doublet
sce <- Seurat::as.SingleCellExperiment(fshd_seurat)
SummarizedExperiment::colData(sce)$sample <- fshd_seurat$sample
set.seed(123)

sce <- scDblFinder::scDblFinder(
  sce,
  samples = "sample"
)
fshd_seurat$scDblFinder.class <- sce$scDblFinder.class
fshd_seurat$scDblFinder.score <- sce$scDblFinder.score
fshd_seurat$scDblFinder.sample <- sce$scDblFinder.sample
fshd_seurat$scDblFinder.weighted <- sce$scDblFinder.weighted
fshd_seurat$scDblFinder.cxds_score <- sce$scDblFinder.cxds_score

fshd_seurat <- subset(
  fshd_seurat,
  subset = scDblFinder.class == "singlet"
)

# QC metrics AFTER filtering
qc_plot_after <- Seurat::VlnPlot(
  fshd_seurat,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "sample",
  ncol = 2,
  pt.size = 0
)


# -------------------------------------------------------------------------------------------

# normalization, highly variable genes, scaling, PCA
fshd_seurat <- Seurat::NormalizeData(fshd_seurat)
fshd_seurat <- Seurat::FindVariableFeatures(
  fshd_seurat,
  selection.method = "vst"
)
fshd_seurat <- Seurat::ScaleData(fshd_seurat)
fshd_seurat <- Seurat::RunPCA(fshd_seurat, npcs = 50)

# inspect PCs before choosing dims: elbow plot
elbow_plot <- SigBridgeR::FindRobustElbow(fshd_seurat, ndims = 50)

# neighborhood graph + clustering (first 16 PCs per the reference workflow)
fshd_seurat <- Seurat::FindNeighbors(fshd_seurat, dims = 1:16)
fshd_seurat <- Seurat::FindClusters(fshd_seurat)

# nonlinear dimensionality reduction for visualization
fshd_seurat <- Seurat::RunTSNE(fshd_seurat, dims = 1:16)
fshd_seurat <- Seurat::RunUMAP(fshd_seurat, dims = 1:16)

# visualization
dim_plot <- Seurat::DimPlot(
  fshd_seurat,
  reduction = "umap",
  group.by = c("sample", "condition", "seurat_clusters"),
  ncol = 3
)


ensembl_id <- rownames(fshd_seurat)

# Genes of bulk data we got were transformed
ens2sym <- IDConverter::convert_hm_genes(ensembl_id, "ensembl")
cli::cli_alert_warning("Found {sum(is.na(ens2sym))} NA value{?s}")

which_na <- which(is.na(ens2sym))
which_dup <- which(duplicated(ens2sym))
ens2sym[union(which_na, which_dup)] <- ensembl_id[
  union(which_na, which_dup)
]

rownames(fshd_seurat) <- ens2sym

qs::qsave(
  fshd_seurat,
  "/home/data/sigbridger/benchmark_data/fshd/fshd_GSE122873_seurat.qs",
  nthreads = 8L
)
