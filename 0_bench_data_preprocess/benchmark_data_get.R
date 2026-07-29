setwd(usethis::proj_path())


library(GEOquery)
library(dplyr)
library(limma)

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
