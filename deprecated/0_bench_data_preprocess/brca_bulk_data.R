#  GSE42568, GSE162228 and TCGA

setwd(usethis::proj_path())
library(rtracklayer)
library(GenomicRanges)

mirai::daemons(8L)

# * log2 GC-RMA signal intensity
# gse42568 <- qs::qread(file.path(data_dir, "brca_bulkdata_GSE42568.qs"))
# gse162228 <- qs::qread(file.path(data_dir, "brca_bulkdata_GSE162228.qs"))

# counts
tcga <- qs::qread(file.path(data_dir, "brca_bulkdata_TCGA.qs"))
source("0_bench_data_preprocess/Counts2TPM.R")
library(TCGAbiolinks)
library(SummarizedExperiment)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
se <- GDCprepare(query)

assayNames(se)
