# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_negative_compare/ov")

set.seed(123)

library(dplyr)
library(data.table)

# ! ov-sc-GSE165897
# ! 阴性对照组

data_path = '/home/data/sigbridger/benchmark_data/ov'
seurat <- qs::qread(
  file.path(data_path, 'hgsoc_GSE165897_seurat.qs'),
  nthreads = 4L
)

# ! bulk
GSE9891 <- qs::qread(
  file.path(data_path, 'ov_bulkdata_GSE9891.qs'),
  nthreads = 4L
)
GSE32062 <- qs::qread(
  file.path(data_path, 'ov_bulkdata_GSE32062_GPL6480.qs'),
  nthreads = 4L
)
GSE140082 <- qs::qread(
  file.path(data_path, 'ov_bulkdata_GSE140082.qs'),
  nthreads = 4L
)


# ! random 20 gene
gene_sc <- rownames(seurat)
gene_b9 <- rownames(GSE9891)
gene_b3 <- rownames(GSE32062)
gene_b12 <- rownames(GSE140082)

common_gene <- Reduce(intersect, list(gene_sc, gene_b9, gene_b3, gene_b12))
message("Total common genes: ", length(common_gene))

# ✅ 重复抽样 100 次，每次 20 个基因
set.seed(123) # 保证每次运行结果一致（若需不同结果可移除此行或换种子）
n_reps <- 100
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

# 转为 data.frame：每列是一次抽样（20 行 × 100 列）
random20_matrix <- do.call(cbind, samples_list)
colnames(random20_matrix) <- paste0("Sample_", 1:n_reps)

# 转为 data.table 并保存
dt_out <- data.table::as.data.table(random20_matrix)
data.table::fwrite(dt_out, file = 'ov_random20_markers_100rep.csv')
