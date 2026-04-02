# ?不考虑
# library(TCGAbiolinks)

# query <- GDCquery(
#     project = "TCGA-LUAD",
#     data.category = "Transcriptome Profiling",
#     data.type = "Gene Expression Quantification"
# )
setwd(usethis::proj_path())


library(GEOquery)
library(dplyr)
library(limma)

luad_bulk_GSE30219 = getGEO("GSE30219")

# filter NA
luad_phenodata_GSE30219 = pData(
  luad_bulk_GSE30219$GSE30219_series_matrix.txt.gz
)

luad_bulkdata_GSE30219 = luad_bulk_GSE30219$GSE30219_series_matrix.txt.gz@assayData$exprs

# translate gene id
luad_feature_GSE30219 = fData(luad_bulk_GSE30219$GSE30219_series_matrix.txt.gz)

gene_symbols <- sapply(
  strsplit(as.character(luad_feature_GSE30219$`Gene Symbol`), " /// "),
  function(x) x[1]
)

rownames(luad_bulkdata_GSE30219) <- gene_symbols

# > any(is.na(luad_bulkdata_GSE30219))
# FALSE

qs::qsave(
  luad_bulkdata_GSE30219,
  "benchmark_data/luad/luad_bulkdata_GSE30219.qs",
  nthreads = 4
)

qs::qsave(
  luad_phenodata_GSE30219,
  "benchmark_data/luad/luad_phenodata_GSE30219.qs",
  nthreads = 4
)

# * peek
# > luad_bulkdata_GSE30219[1:4,1:4]
#           GSM748053 GSM748054 GSM748055 GSM748056
# 1007_s_at 10.805665 10.381745 11.119014 10.169954
# 1053_at    6.910219  7.921614  8.172580  7.364354
# 117_at     6.058989  5.964559  5.944449  6.061103
# 121_at     7.147779  7.416915  7.683858  7.206680

# > dim(luad_bulkdata_GSE30219)
# [1] 54675   307

ad_bulk_GSE39420 = getGEO("GSE39420")

ad_phenodata_GSE39420 = pData(ad_bulk_GSE39420$GSE39420_series_matrix.txt.gz)

ad_bulkdata_GSE39420 = ad_bulk_GSE39420$GSE39420_series_matrix.txt.gz@assayData$exprs

ad_feature_GSE39420 = fData(ad_bulk_GSE39420$GSE39420_series_matrix.txt.gz)

# for gene name convertion
# platform <- annotation(ad_bulk_GSE39420[[1]])

gene_symbols <- ifelse(
  grepl(" /// ", ad_feature_GSE39420[["gene_assignment"]], fixed = TRUE),
  sapply(
    strsplit(
      ad_feature_GSE39420[["gene_assignment"]],
      " /// ",
      fixed = TRUE
    ),
    `[`,
    1
  ),
  rownames(ad_feature_GSE39420)
)

rownames(ad_bulkdata_GSE39420) <- gene_symbols

qs::qsave(
  ad_bulkdata_GSE39420,
  "benchmark_data/ad/ad_bulkdata_GSE39420.qs",
  nthreads = 4
)

qs::qsave(
  ad_phenodata_GSE39420,
  "benchmark_data/ad/ad_phenodata_GSE39420.qs",
  nthreads = 4
)

# > ad_bulkdata_GSE39420[1:4,1:4]
#         GSM967918 GSM967919 GSM967920 GSM967921
# 7892501     2.565     2.703     2.425     2.222
# 7892502     3.379     3.780     3.660     3.437
# 7892503     2.429     2.752     2.441     2.780
# 7892504     9.329     8.628     8.017     8.522

ad_bulk_GSE109887 = getGEO("GSE109887")

ad_bulkdata_GSE109887 = ad_bulk_GSE109887$GSE109887_series_matrix.txt.gz@assayData$exprs

ad_pheno_GSE109887 = pData(ad_bulk_GSE109887$GSE109887_series_matrix.txt.gz)

# * No need to translate gene names, already converted
# ad_feature_GSE109887 = fData(ad_bulk_GSE109887$GSE109887_series_matrix.txt.gz)

qs::qsave(
  ad_bulkdata_GSE109887,
  "benchmark_data/ad/ad_bulkdata_GSE109887.qs",
  nthreads = 4
)
qs::qsave(
  ad_pheno_GSE109887,
  "benchmark_data/ad/ad_pheno_GSE109887.qs",
  nthreads = 4
)


ad_bulk_GSE28146 = getGEO("GSE28146")

ad_bulkdata_GSE28146 = ad_bulk_GSE28146$GSE28146_series_matrix.txt.gz@assayData$exprs

ad_pheno_GSE28146 = pData(ad_bulk_GSE28146$GSE28146_series_matrix.txt.gz)

ad_feature_GSE28146 = fData(ad_bulk_GSE28146$GSE28146_series_matrix.txt.gz)

# convert gene symbols
gene_symbols <- sapply(
  strsplit(as.character(ad_feature_GSE28146$`Gene Symbol`), " /// "),
  function(x) x[1]
)
rownames(ad_bulkdata_GSE28146) <- gene_symbols

qs::qsave(
  ad_bulkdata_GSE28146,
  "benchmark_data/ad/ad_bulkdata_GSE28146.qs",
  nthreads = 4
)
qs::qsave(
  ad_pheno_GSE28146,
  "benchmark_data/ad/ad_pheno_GSE28146.qs",
  nthreads = 4
)

# ! It is empty via downloading from online or GEO dataset
# ! OK
# brca_sc_GSE161529 = getGEO(
#     "GSE161529"
# )

# getGEOSuppFiles("GSE161529")

# brca_sc_mat_GSE161529 = brca_sc_GSE161529$GSE161529_series_matrix.txt.gz@assayData$exprs
# sc = GSMList(brca_sc_GSE161529)

brca_sc_mat_GSE161529 = data.table::fread(
  "benchmark_data/brca/GSE161529_features.tsv"
)

# -----------------------------
brca_bulk_GSE162228 = getGEO("GSE162228")

brca_bulkdata_GSE162228 = brca_bulk_GSE162228$GSE162228_series_matrix.txt.gz@assayData$exprs

brca_pheno_GSE162228 = pData(brca_bulk_GSE162228$GSE162228_series_matrix.txt.gz)

brca_feature_GSE162228 = fData(
  brca_bulk_GSE162228$GSE162228_series_matrix.txt.gz
)
#
gene_symbols <- sapply(
  strsplit(as.character(brca_feature_GSE162228$`Gene Symbol`), " /// "),
  function(x) x[1]
)
rownames(brca_bulkdata_GSE162228) <- gene_symbols

qs::qsave(
  brca_bulkdata_GSE162228,
  "benchmark_data/brca/brca_bulkdata_GSE162228.qs",
  nthreads = 4
)
qs::qsave(
  brca_pheno_GSE162228,
  "benchmark_data/brca/brca_pheno_GSE162228.qs",
  nthreads = 4
)


brca_bulk_GSE42568 = getGEO("GSE42568")

brca_bulkdata_GSE42568 = brca_bulk_GSE42568$GSE42568_series_matrix.txt.gz@assayData$exprs

brca_pheno_GSE42568 = pData(brca_bulk_GSE42568$GSE42568_series_matrix.txt.gz)

brca_feature_GSE42568 = fData(brca_bulk_GSE42568$GSE42568_series_matrix.txt.gz)

gene_symbols <- sapply(
  strsplit(as.character(brca_feature_GSE42568$`Gene Symbol`), " /// "),
  function(x) x[1]
)
rownames(brca_bulkdata_GSE42568) <- gene_symbols

qs::qsave(
  brca_bulkdata_GSE42568,
  "benchmark_data/brca/brca_bulkdata_GSE42568.qs",
  nthreads = 4
)
qs::qsave(
  brca_pheno_GSE42568,
  "benchmark_data/brca/brca_pheno_GSE42568.qs",
  nthreads = 4
)
# -----------------------------
# ! it is empty via downloading from online

# ov_sc_GSE184880 = getGEO("GSE184880")

# ov_scmat_GSE184880 = ov_sc_GSE184880$GSE184880_series_matrix.txt.gz@assayData$exprs

# ov_pheno_GSE184880 = pData(ov_sc_GSE184880$GSE184880_series_matrix.txt.gz)

# ov_feature_GSE184880 = fData(ov_sc_GSE184880$GSE184880_series_matrix.txt.gz)

# gene_symbols <- sapply(
#     strsplit(as.character(ov_feature_GSE184880$`Gene Symbol`), " /// "),
#     function(x) x[1]
# )
# rownames(ov_scmat_GSE184880) <- gene_symbols

# qs::qsave(
#     ov_scmat_GSE184880,
#     "benchmark_data/ov/ov_scmat_GSE184880.qs",
#     nthreads = 4
# )
# qs::qsave(
#     ov_pheno_GSE184880,
#     "benchmark_data/ov/ov_pheno_GSE184880.qs",
#     nthreads = 4
# )

# ! it is empty via downloading from online
# ! OK

# ov_sc_GSE165897 = getGEO("GSE165897")

# ov_scmat_GSE165897 = exprs(
#     ov_sc_GSE165897$`GSE165897-GPL16791_series_matrix.txt.gz`
# )

# * unnecessary, it's the single-cell raw matrix
# ov_sc_GSE165897 = data.table::fread(
#     "benchmark_data/ov/GSE165897_UMIcounts_HGSOC.tsv"
# ) #<<<HERE

# ov_sc_GSE165897 = ov_sc_GSE165897 %>% tibble::column_to_rownames("V1")

# ---------------------------

ov_bulk_GSE9891 = getGEO("GSE9891")

ov_bulkdata_GSE9891 = exprs(ov_bulk_GSE9891$GSE9891_series_matrix.txt.gz)
ov_pheno_GSE9891 = pData(ov_bulk_GSE9891$GSE9891_series_matrix.txt.gz)
ov_feature_GSE9891 = fData(ov_bulk_GSE9891$GSE9891_series_matrix.txt.gz)

gene_symbols <- sapply(
  strsplit(as.character(ov_feature_GSE9891$`Gene Symbol`), " /// "),
  function(x) x[1]
)
rownames(ov_bulkdata_GSE9891) <- gene_symbols

qs::qsave(
  ov_bulkdata_GSE9891,
  "benchmark_data/ov/ov_bulkdata_GSE9891.qs",
  nthreads = 4
)
qs::qsave(
  ov_pheno_GSE9891,
  "benchmark_data/ov/ov_pheno_GSE9891.qs",
  nthreads = 4
)

ov_bulk_GSE140082 = getGEO("GSE140082")

ov_bulkdata_GSE140082 = exprs(ov_bulk_GSE140082$GSE140082_series_matrix.txt.gz)
ov_pheno_GSE140082 = pData(ov_bulk_GSE140082$GSE140082_series_matrix.txt.gz)
ov_feature_GSE140082 = fData(ov_bulk_GSE140082$GSE140082_series_matrix.txt.gz)

gene_symbols <- sapply(
  strsplit(as.character(ov_feature_GSE140082$`Symbol`), " /// "),
  function(x) x[1]
)
rownames(ov_bulkdata_GSE140082) <- gene_symbols

qs::qsave(
  ov_bulkdata_GSE140082,
  "benchmark_data/ov/ov_bulkdata_GSE140082.qs",
  nthreads = 4
)
qs::qsave(
  ov_pheno_GSE140082,
  "benchmark_data/ov/ov_pheno_GSE140082.qs",
  nthreads = 4
)


# GSE32062
# ! GPL570 was not used, GPL6480 was used
# ! OK
ov_bulk_GSE32062 = getGEO("GSE32062")

ov_bulkdata_GSE32062 = exprs(
  ov_bulk_GSE32062$`GSE32062-GPL570_series_matrix.txt.gz`
)
ov_pheno_GSE32062 = pData(
  ov_bulk_GSE32062$`GSE32062-GPL570_series_matrix.txt.gz`
)
ov_feature_GSE32062 = fData(
  ov_bulk_GSE32062$`GSE32062-GPL570_series_matrix.txt.gz`
)

gene_symbols <- sapply(
  strsplit(as.character(ov_feature_GSE32062$`Gene Symbol`), " /// "),
  function(x) x[1]
)
rownames(ov_bulkdata_GSE32062) <- gene_symbols

qs::qsave(
  ov_bulkdata_GSE32062,
  "benchmark_data/ov/ov_bulkdata_GSE32062_GPL570.qs",
  nthreads = 4
)
qs::qsave(
  ov_pheno_GSE32062,
  "benchmark_data/ov/ov_pheno_GSE32062_GPL570.qs",
  nthreads = 4
)

# ! partial empty data
# ! OK
ov_bulk_GSE32062$`GSE32062-GPL6480_series_matrix.txt.gz` |> fData() |> View()
ov_bulk_GSE32062$`GSE32062-GPL6480_series_matrix.txt.gz` |> pData() |> View()

# ! empty via `getGEO`
# ! OK, but duplicated genes exist
# kidney_sc_GSE175540 = getGEO("GSE175540")

# kidney_scmat_GSE175540 = exprs(
#     kidney_sc_GSE175540$GSE175540_series_matrix.txt.gz
# )
kidney_sc_GSE175540 = data.table::fread(
  "benchmark_data/kidney/GSE175540_batch_corrected_bulk_RNA_seq_TPM_25_02_22.csv"
)
kidney_sc_GSE175540 = kidney_sc_GSE175540 %>% tibble::column_to_rownames("V1")


# ! empty via `getGEO`
colon_sc_GSE225857 = getGEO("GSE225857")

colon_scmat_GSE225857 = exprs(
  colon_sc_GSE225857$GSE225857_series_matrix.txt.gz
)
