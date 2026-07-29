#  GSE42568, GSE162228 and TCGA

setwd(usethis::proj_path())
library(rtracklayer)
library(GenomicRanges)

mirai::daemons(8L)

# * log2 GC-RMA signal intensity
# gse42568 <- qs::qread(file.path(data_dir, "brca_bulkdata_GSE42568.qs"))
# gse162228 <- qs::qread(file.path(data_dir, "brca_bulkdata_GSE162228.qs"))

# counts
data_dir <- "/home/data/sigbridger/benchmark_data/brca"
tcga <- qs::qread(file.path(data_dir, "brca_bulkdata_TCGA.qs"))
source("0_bench_data_preprocess/Counts2TPM.R")

gtf_file <- "/home/data/data-resource/TCGA/gencode.v33.annotation.gtf.gz"
gtf <- rtracklayer::import(gtf_file)

# 只取 exon
exons <- gtf[gtf$type == "exon"]

# 按 gene_id 分组
exons_by_gene <- split(exons, exons$gene_id)

# 对每个 gene 的 exon 做 reduce，去掉重叠部分，然后求长度
reduced_exons_by_gene <- purrr::map(
  exons_by_gene,
  purrr::in_parallel(\(x) IRanges::reduce(x)),
  .progress = "Reducing exons by gene"
)
gene_length <- purrr::map(
  reduced_exons_by_gene,
  purrr::in_parallel(\(x) sum(IRanges::width(x))),
  .progress = "Calculating gene length"
)

gene_length_df <- data.frame(
  gene_id = names(gene_length),
  gene_length = as.numeric(gene_length)
)


load("/home/data/data-resource/TCGA/gencode.v33.annotation.Rdata") # gencode.v33.annotation
gencode.v33.annotation <- gencode.v33.annotation[
  !is.na(gencode.v33.annotation$Ensembl),
]
if (all(gencode.v33.annotation$Ensembl == gene_length_df$gene_id)) {
  print("gene_id match")
}

gene_length_cmb <- cbind(gencode.v33.annotation, gene_length_df)

data.table::fwrite(
  gene_length_cmb,
  "0_bench_data_preprocess/gencode_v33_tcga_gene_length.csv"
)

gene_length_list <- setNames(
  gene_length_cmb$gene_length,
  gene_length_cmb$gene_id
)
tcga_tpm <- Counts2TPM(as.matrix(tcga), gene_length_list)
qs::qsave(
  tcga_tpm,
  file.path(data_dir, "brca_bulkdata_TCGA_tpm.qs"),
  nthreads = 4L
)

mirai::daemons(0L)
gc()
