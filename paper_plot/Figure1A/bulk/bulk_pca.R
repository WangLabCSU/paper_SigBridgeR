setwd(file.path(
  usethis::proj_path(),
  "paper_plot/Figure1A/bulk"
))

source("DrawBulkPCA.R")

data_dir <- "/home/data/sigbridger/benchmark_data/"

pheno <- qs::qread(file.path(data_dir, "lung/TCGA_LUAD_pheno.qs"))
tpm <- qs::qread(file.path(data_dir, "lung/TCGA_LUAD_bulkdata_tpm.qs"))

cm_samples <- intersect(rownames(pheno), colnames(tpm))

tpm <- tpm[, cm_samples]
pheno <- pheno[cm_samples, ]


# log 化压缩动态范围
tpm_log <- log2(tpm + 1)

# 特征选择：取方差最大的 2000 个基因。全部 5.7 万个基因里有 4 个零方差基因，
# scale. = TRUE 会因此产生 NaN。
gene_var <- matrixStats::rowVars(tpm_log)
tpm_hvg <- tpm_log[order(gene_var, decreasing = TRUE)[seq_len(2000L)], ]

pca <- DrawBulkPCA(tpm_hvg, group = pheno$tissue_type)

# 移除 DrawBulkPCA() 添加的 patchwork 总标题。
pca$patches$annotation$title <- NULL

pca

ggplot2::ggsave("bulk_pca.png", pca, width = 10, height = 10)
