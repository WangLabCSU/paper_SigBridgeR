setwd(usethis::proj_path())


library(GEOquery)
library(dplyr)
library(limma)

source(
  "0_bench_data_preprocess/convert_gene_symbol.R"
)
source(
  "0_bench_data_preprocess/fpkm_to_tpm.R"
)

out_dir <- "/home/data/sigbridger/"

# -----------------------------
brca_bulk_GSE162228 <- GEOquery::getGEO("GSE162228")

brca_bulkdata_GSE162228 <- brca_bulk_GSE162228$GSE162228_series_matrix.txt.gz@assayData$exprs

brca_pheno_GSE162228 <- pData(
  brca_bulk_GSE162228$GSE162228_series_matrix.txt.gz
)

brca_feature_GSE162228 <- fData(
  brca_bulk_GSE162228$GSE162228_series_matrix.txt.gz
)

gene_symbols <- sapply(
  strsplit(as.character(brca_feature_GSE162228$`Gene Symbol`), " /// "),
  function(x) x[1]
)
rownames(brca_bulkdata_GSE162228) <- gene_symbols

qs::qsave(
  brca_bulkdata_GSE162228,
  file.path(out_dir, "benchmark_data/brca/brca_bulkdata_GSE162228.qs"),
  nthreads = 4
)
qs::qsave(
  brca_pheno_GSE162228,
  file.path(out_dir, "benchmark_data/brca/brca_pheno_GSE162228.qs"),
  nthreads = 4
)


brca_bulk_GSE42568 <- GEOquery::getGEO("GSE42568")

brca_bulkdata_GSE42568 <- brca_bulk_GSE42568$GSE42568_series_matrix.txt.gz@assayData$exprs

brca_pheno_GSE42568 <- pData(brca_bulk_GSE42568$GSE42568_series_matrix.txt.gz)

brca_feature_GSE42568 <- fData(brca_bulk_GSE42568$GSE42568_series_matrix.txt.gz)

gene_symbols <- sapply(
  strsplit(as.character(brca_feature_GSE42568$`Gene Symbol`), " /// "),
  function(x) x[1]
)
rownames(brca_bulkdata_GSE42568) <- gene_symbols

qs::qsave(
  brca_bulkdata_GSE42568,
  file.path(out_dir, "benchmark_data/brca/brca_bulkdata_GSE42568.qs"),
  nthreads = 4
)
qs::qsave(
  brca_pheno_GSE42568,
  file.path(out_dir, "benchmark_data/brca/brca_pheno_GSE42568.qs"),
  nthreads = 4
)


# ---------------------------

ov_bulk_GSE9891 <- GEOquery::getGEO("GSE9891")

ov_bulkdata_GSE9891 <- exprs(ov_bulk_GSE9891$GSE9891_series_matrix.txt.gz)
ov_pheno_GSE9891 <- pData(ov_bulk_GSE9891$GSE9891_series_matrix.txt.gz)
ov_feature_GSE9891 <- fData(ov_bulk_GSE9891$GSE9891_series_matrix.txt.gz)

gene_symbols <- sapply(
  strsplit(as.character(ov_feature_GSE9891$`Gene Symbol`), " /// "),
  function(x) x[1]
)
rownames(ov_bulkdata_GSE9891) <- gene_symbols

qs::qsave(
  ov_bulkdata_GSE9891,
  file.path(out_dir, "benchmark_data/ov/ov_bulkdata_GSE9891.qs"),
  nthreads = 4
)
qs::qsave(
  ov_pheno_GSE9891,
  file.path(out_dir, "benchmark_data/ov/ov_pheno_GSE9891.qs"),
  nthreads = 4
)

ov_bulk_GSE140082 <- GEOquery::getGEO("GSE140082")

ov_bulkdata_GSE140082 <- exprs(ov_bulk_GSE140082$GSE140082_series_matrix.txt.gz)
ov_pheno_GSE140082 <- pData(ov_bulk_GSE140082$GSE140082_series_matrix.txt.gz)
ov_feature_GSE140082 <- fData(ov_bulk_GSE140082$GSE140082_series_matrix.txt.gz)

gene_symbols <- sapply(
  strsplit(as.character(ov_feature_GSE140082$`Symbol`), " /// "),
  function(x) x[1]
)
rownames(ov_bulkdata_GSE140082) <- gene_symbols

qs::qsave(
  ov_bulkdata_GSE140082,
  file.path(out_dir, "benchmark_data/ov/ov_bulkdata_GSE140082.qs"),
  nthreads = 4
)
qs::qsave(
  ov_pheno_GSE140082,
  file.path(out_dir, "benchmark_data/ov/ov_pheno_GSE140082.qs"),
  nthreads = 4
)


# GSE32062
# ! GPL570 was not used, GPL6480 was used
# ! OK
ov_bulk_GSE32062 <- GEOquery::getGEO("GSE32062")

# ov_bulkdata_GSE32062 <- exprs(
#   ov_bulk_GSE32062$`GSE32062-GPL570_series_matrix.txt.gz`
# )
# ov_pheno_GSE32062 <- pData(
#   ov_bulk_GSE32062$`GSE32062-GPL570_series_matrix.txt.gz`
# )
# ov_feature_GSE32062 <- fData(
#   ov_bulk_GSE32062$`GSE32062-GPL570_series_matrix.txt.gz`
# )

# gene_symbols <- sapply(
#   strsplit(as.character(ov_feature_GSE32062$`Gene Symbol`), " /// "),
#   function(x) x[1]
# )
# rownames(ov_bulkdata_GSE32062) <- gene_symbols

# qs::qsave(
#   ov_bulkdata_GSE32062,
#   file.path(out_dir, "benchmark_data/ov/ov_bulkdata_GSE32062_GPL570.qs"),
#   nthreads = 4
# )
# qs::qsave(
#   ov_pheno_GSE32062,
#   file.path(out_dir, "benchmark_data/ov/ov_pheno_GSE32062_GPL570.qs"),
#   nthreads = 4
# )

ov_bulkdata_GSE32062 <- exprs(
  ov_bulk_GSE32062$`GSE32062-GPL6480_series_matrix.txt.gz`
)

ov_pheno_GSE32062 <- pData(
  ov_bulk_GSE32062$`GSE32062-GPL6480_series_matrix.txt.gz`
)
ov_feature_GSE32062 <- fData(
  ov_bulk_GSE32062$`GSE32062-GPL6480_series_matrix.txt.gz`
)

gene_symbols <- sapply(
  strsplit(as.character(ov_feature_GSE32062$GENE_SYMBOL), " /// "),
  function(x) x[1]
)
# gene must match gene symbols in single cell data, so probes with no gene symbol are removed
ov_bulkdata_GSE32062 <- ov_bulkdata_GSE32062[!is.na(gene_symbols), ]
rownames(ov_bulkdata_GSE32062) <- gene_symbols[!is.na(gene_symbols)]

qs::qsave(
  ov_bulkdata_GSE32062,
  file.path(out_dir, "benchmark_data/ov/ov_bulkdata_GSE32062_GPL6480.qs"),
  nthreads = 4
)
qs::qsave(
  ov_pheno_GSE32062,
  file.path(out_dir, "benchmark_data/ov/ov_pheno_GSE32062.qs"),
  nthreads = 4
)
# -------------------------------------------------------------------------------------------

ad_bulk_GSE39420 <- geokit::geo_matrix(
  "GSE39420",
  odir = "/home/data/sigbridger/benchmark_data/ad/GSE39420",
  add_gpl = TRUE
)
ad_bulkdata_GSE39420 <- Biobase::exprs(ad_bulk_GSE39420)
ad_pheno_GSE39420 <- Biobase::pData(ad_bulk_GSE39420)
ad_feature_GSE39420 <- Biobase::fData(ad_bulk_GSE39420)

which_na <- is.na(ad_feature_GSE39420$`Gene symbol`)

ad_bulkdata_GSE39420 <- ad_bulkdata_GSE39420[!which_na, ]
ad_feature_GSE39420 <- ad_feature_GSE39420[!which_na, ]
rownames(ad_bulkdata_GSE39420) <- ad_feature_GSE39420$`Gene symbol` |>
  stringr::str_remove("//.*$") |>
  stringr::str_trim()

# * kept controls and early-onset Alzheimer’s disease patients to identify the corresponding cell subpopulations
which_ctrl_earlyonset <- stringr::str_detect(
  ad_pheno_GSE39420$title,
  "control|onset"
)
ad_bulkdata_GSE39420 <- ad_bulkdata_GSE39420[, which_ctrl_earlyonset]
ad_pheno_GSE39420 <- ad_pheno_GSE39420[which_ctrl_earlyonset, ]

qs::qsave(
  ad_bulkdata_GSE39420,
  file.path(out_dir, "benchmark_data/ad/ad_bulkdata_GSE39420.qs"),
  nthreads = 4
)
qs::qsave(
  ad_pheno_GSE39420,
  file.path(out_dir, "benchmark_data/ad/ad_pheno_GSE39420.qs"),
  nthreads = 4
)

# --------------------------------------------------------------------------------------------------------
fshd_bulk_GSE140261 <- geokit::geo_matrix(
  "GSE140261",
  odir = "/home/data/sigbridger/benchmark_data/fshd/GSE140261",
  add_gpl = TRUE
)
fshd_supple_GSE140261 <- geokit::geo_suppl(
  "GSE140261",
  odir = "/home/data/sigbridger/benchmark_data/fshd/GSE140261"
)
# fshd_bulkdata_GSE140261 <- Biobase::exprs(fshd_bulk_GSE140261) # NULL
fshd_pheno_GSE140261 <- Biobase::pData(fshd_bulk_GSE140261)
# fshd_feature_GSE140261 <- Biobase::fData(fshd_bulk_GSE140261) # NULL

fshd_bulkdata_GSE140261 <- data.table::fread(
  "/home/data/sigbridger/benchmark_data/fshd/GSE140261/GSE140261_year2_TPM.csv"
)
genes <- fshd_bulkdata_GSE140261[, 1][[1]]
gene_symbols <- convert_gene_symbol(genes)

fshd_bulkdata_GSE140261 <- as.matrix(fshd_bulkdata_GSE140261[, -1])
rownames(fshd_bulkdata_GSE140261) <- gene_symbols

qs::qsave(
  fshd_bulkdata_GSE140261,
  file.path(out_dir, "benchmark_data/fshd/fshd_bulkdata_GSE140261.qs"),
  nthreads = 4
)
qs::qsave(
  fshd_pheno_GSE140261,
  file.path(out_dir, "benchmark_data/fshd/fshd_pheno_GSE140261.qs"),
  nthreads = 4
)

# --------------------------------------------------------------------------------------------------------

fshd_sc_GSE122873 <- geokit::geo_matrix(
  "GSE122873",
  odir = "/home/data/sigbridger/benchmark_data/fshd/GSE122873",
  add_gpl = TRUE
)
# fshd_scdata_GSE122873 <- Biobase::exprs(fshd_sc_GSE122873) # NULL
fshd_meta_GSE122873 <- Biobase::pData(fshd_sc_GSE122873)
# fshd_feature_GSE122873 <- Biobase::fData(fshd_sc_GSE122873) # NULL

fshd_supple_GSE122873 <- geokit::geo_suppl(
  "GSE122873",
  odir = "/home/data/sigbridger/benchmark_data/fshd/GSE122873"
)

# --------------------------------------------------------------------------------------------------------
# FSHD validation datasets (bulk)

# * GSE115650
fshd_bulk_GSE115650 <- geokit::geo_matrix(
  "GSE115650",
  odir = "/home/data/sigbridger/benchmark_data/fshd/GSE115650",
  add_gpl = TRUE
)
fshd_supple_GSE115650 <- geokit::geo_suppl(
  "GSE115650",
  odir = "/home/data/sigbridger/benchmark_data/fshd/GSE115650"
)
# fshd_bulkdata_GSE115650 <- Biobase::exprs(fshd_bulk_GSE115650) # NULL
fshd_pheno_GSE115650 <- Biobase::pData(fshd_bulk_GSE115650)
# fshd_feature_GSE115650 <- Biobase::fData(fshd_bulk_GSE115650) # NULL
fshd_bulkdata_GSE115650 <- data.table::fread(
  "/home/data/sigbridger/benchmark_data/fshd/GSE115650/GSE115650_biopsy_FPKM.csv"
)

genes <- fshd_bulkdata_GSE115650[, 1][[1]]
gene_symbols <- convert_gene_symbol(genes)

fshd_bulkdata_GSE115650 <- as.matrix(fshd_bulkdata_GSE115650[, -1])
rownames(fshd_bulkdata_GSE115650) <- gene_symbols

fshd_bulkdata_GSE115650 <- fpkm_to_tpm(fshd_bulkdata_GSE115650)

qs::qsave(
  fshd_bulkdata_GSE115650,
  file.path(out_dir, "benchmark_data/fshd/fshd_bulkdata_GSE115650.qs"),
  nthreads = 4
)
qs::qsave(
  fshd_pheno_GSE115650,
  file.path(out_dir, "benchmark_data/fshd/fshd_pheno_GSE115650.qs"),
  nthreads = 4
)

# ------------------------------------------------------------------------------------------------------

# * GSE56787
fshd_bulk_GSE56787 <- geokit::geo_matrix(
  "GSE56787",
  odir = "/home/data/sigbridger/benchmark_data/fshd/GSE56787",
  add_gpl = TRUE
)
fshd_supple_GSE56787 <- geokit::geo_suppl(
  "GSE56787",
  odir = "/home/data/sigbridger/benchmark_data/fshd/GSE56787"
)
# fshd_bulkdata_GSE56787 <- Biobase::exprs(fshd_bulk_GSE56787) # NULL
fshd_pheno_GSE56787 <- Biobase::pData(fshd_bulk_GSE56787)
# fshd_feature_GSE56787 <- Biobase::fData(fshd_bulk_GSE56787) # NULL

fshd_bulkdata_GSE56787 <- data.table::fread(
  "/home/data/sigbridger/benchmark_data/fshd/GSE56787/GSE56787_biopsy.gene.counts.csv"
)
gene <- fshd_bulkdata_GSE56787[, 1] %>%
  unlist() %>%
  stringr::str_remove_all("\\..*$")
fshd_bulkdata_GSE56787 <- fshd_bulkdata_GSE56787[, -1]
fshd_bulkdata_GSE56787 <- as.matrix(fshd_bulkdata_GSE56787)
rownames(fshd_bulkdata_GSE56787) <- gene

Rcpp::sourceCpp("0_bench_data_preprocess/gtf_lines_to_gene_length.cpp")
gene_length <- gtf_file_to_gene_length(
  "/home/data/sigbridger/benchmark_data/Homo_sapiens.GRCh37.75.gtf.gz"
)

source("0_bench_data_preprocess/counts_to_tpm.R")
fshd_bulkdata_GSE56787 <- counts_to_tpm(fshd_bulkdata_GSE56787, gene_length)

source("0_bench_data_preprocess/convert_gene_symbol.R")
gene_symbols <- convert_gene_symbol(rownames(fshd_bulkdata_GSE56787))
rownames(fshd_bulkdata_GSE56787) <- gene_symbols

# * use the muscle biopsies
muscle_sample <- stringr::str_detect(
  fshd_pheno_GSE56787$title,
  "[A-Z]{1}[0-9]{1,2}$"
)
fshd_pheno_GSE56787 <- fshd_pheno_GSE56787[muscle_sample, ]
# * order
sample_map <- setNames(rownames(fshd_pheno_GSE56787), fshd_pheno_GSE56787$title)
colnames(fshd_bulkdata_GSE56787) <- sample_map[colnames(fshd_bulkdata_GSE56787)]

qs::qsave(
  fshd_bulkdata_GSE56787,
  file.path(out_dir, "benchmark_data/fshd/fshd_bulkdata_GSE56787.qs"),
  nthreads = 4
)
qs::qsave(
  fshd_pheno_GSE56787,
  file.path(out_dir, "benchmark_data/fshd/fshd_pheno_GSE56787.qs"),
  nthreads = 4
)


# ------------------------------------------------------------------------------------------------------
# AD validation datasets (bulk)

# * GSE109887
ad_bulk_GSE109887 <- geokit::geo_matrix(
  "GSE109887",
  odir = "/home/data/sigbridger/benchmark_data/ad/GSE109887",
  add_gpl = TRUE
)
ad_supple_GSE109887 <- geokit::geo_suppl(
  "GSE109887",
  odir = "/home/data/sigbridger/benchmark_data/ad/GSE109887"
)

ad_bulkdata_GSE109887 <- Biobase::exprs(ad_bulk_GSE109887) # already normalized and labelled
ad_pheno_GSE109887 <- Biobase::pData(ad_bulk_GSE109887)
# ad_feature_GSE109887 <- Biobase::fData(ad_bulk_GSE109887) # not used

qs::qsave(
  ad_bulkdata_GSE109887,
  file.path(out_dir, "benchmark_data/ad/ad_bulkdata_GSE109887.qs"),
  nthreads = 4
)
qs::qsave(
  ad_pheno_GSE109887,
  file.path(out_dir, "benchmark_data/ad/ad_pheno_GSE109887.qs"),
  nthreads = 4
)

# ------------------------------------------------------------------------------------------------------

# * GSE28146
ad_bulk_GSE28146 <- geokit::geo_matrix(
  "GSE28146",
  odir = "/home/data/sigbridger/benchmark_data/ad/GSE28146",
  add_gpl = TRUE
)
ad_supple_GSE28146 <- geokit::geo_suppl(
  "GSE28146",
  odir = "/home/data/sigbridger/benchmark_data/ad/GSE28146"
)

ad_bulkdata_GSE28146 <- Biobase::exprs(ad_bulk_GSE28146)
ad_pheno_GSE28146 <- Biobase::pData(ad_bulk_GSE28146)
ad_feature_GSE28146 <- Biobase::fData(ad_bulk_GSE28146)

gene_symbols <- stringr::str_remove(ad_feature_GSE28146$`Gene symbol`, "///.*$")
# anyDuplicated(gene_symbols)
ad_bulkdata_GSE28146 <- aggregate(
  ad_bulkdata_GSE28146,
  by = list(gene = gene_symbols),
  FUN = sum
)
ad_bulkdata_GSE28146 <- tibble::column_to_rownames(ad_bulkdata_GSE28146, "gene")
ad_bulkdata_GSE28146 <- log2(ad_bulkdata_GSE28146)

qs::qsave(
  ad_bulkdata_GSE28146,
  file.path(out_dir, "benchmark_data/ad/ad_bulkdata_GSE28146.qs"),
  nthreads = 4
)
qs::qsave(
  ad_pheno_GSE28146,
  file.path(out_dir, "benchmark_data/ad/ad_pheno_GSE28146.qs"),
  nthreads = 4
)
