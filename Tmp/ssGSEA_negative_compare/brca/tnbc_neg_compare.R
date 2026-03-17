# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_negative_compare/brca")

set.seed(123)

library(dplyr)
library(data.table)

# ! luad-sc-161529-tnbc
# ! 阳性对照组

data_path = '/home/data/sigbridger/benchmark_data/brca'
seurat <- qs::qread(
  file.path(data_path, 'seurat_tnbc.qs'),
  nthreads = 4
)


# ! bulk
GSE42568 <- qs::qread(
  file.path(data_path, 'brca_bulkdata_GSE42568.qs'),
  nthreads = 4
)
GSE162228 <- qs::qread(
  file.path(data_path, 'brca_bulkdata_GSE162228.qs'),
  nthreads = 4
)
TCGA_BRCA <- qs::qread(
  file.path(data_path, 'brca_bulkdata_TCGA.qs'),
  nthreads = 4
)


# ! random 20 gene
gene_sc <- rownames(seurat)
gene_b4 <- rownames(GSE42568)
gene_b1 <- rownames(GSE162228)
gene_bt <- rownames(TCGA_BRCA)

common_gene <- Reduce(
  intersect,
  list(gene_sc, gene_b4, gene_b1, gene_bt)
)
message("Total common genes: ", length(common_gene))

# ✅ 重复抽样 10 次，每次 20 个基因
set.seed(123) # 保证每次运行结果一致（若需不同结果可移除此行或换种子）
n_reps <- 10
n_genes <- 20

# 检查：确保 common_gene 足够多
if (length(common_gene) < n_genes) {
  stop("Not enough common genes to sample 20 unique genes!")
}

# 方法1：用 replicate + simplify = FALSE → 转为 data.frame 列
samples_list <- replicate(
  n = n_reps,
  expr = sample(common_gene, size = n_genes, replace = FALSE),
  simplify = FALSE
)

# 转为 data.frame：每列是一次抽样（20 行 × 10 列）
random20_matrix <- do.call(cbind, samples_list)
colnames(random20_matrix) <- paste0("Sample_", 1:n_reps)

# 转为 data.table 并保存
dt_out <- data.table::as.data.table(random20_matrix)
data.table::fwrite(dt_out, file = 'tnbc_random20_markers_10rep.csv')
