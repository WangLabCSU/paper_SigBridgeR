# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(file.path(usethis::proj_path(), "3_negative_ctrl/luad"))

set.seed(123)

library(dplyr)
library(data.table)

# ! ov-sc-GSE165897

data_path = '/home/data/sigbridger/benchmark_data/lung'
seurat <- qs::qread(
  file.path(data_path, 'luad_GSE123902_seurat.qs'),
  nthreads = 4
)

# ! bulk
GSE3141 <- qs::qread(
  file.path(data_path, 'lung_bulkdata_GSE3141.qs'),
  nthreads = 4
)
GSE8894 <- qs::qread(
  file.path(data_path, 'lung_bulkdata_GSE8894.qs'),
  nthreads = 4
)
GSE31210 <- qs::qread(
  file.path(data_path, 'lung_bulkdata_GSE31210.qs'),
  nthreads = 4
)
TCGA_LUAD <- qs::qread(
  file.path(data_path, 'TCGA_LUAD_bulkdata.qs'),
  nthreads = 4
)


# ! random 20 gene
gene_sc <- rownames(seurat)
gene_b3 <- rownames(GSE3141)
gene_b8 <- rownames(GSE8894)
gene_31 <- rownames(GSE31210)
gene_bt <- rownames(TCGA_LUAD)

common_gene <- Reduce(
  intersect,
  list(gene_sc, gene_b3, gene_b8, gene_31, gene_bt)
)
message("Total common genes: ", length(common_gene))


set.seed(123)
n_reps <- 100
n_genes <- 20

if (length(common_gene) < n_genes) {
  stop("Not enough common genes to sample 20 unique genes!")
}

samples_list <- replicate(
  n = n_reps,
  expr = sample(common_gene, size = n_genes, replace = FALSE),
  simplify = FALSE
)

random20_matrix <- do.call(cbind, samples_list)
colnames(random20_matrix) <- paste0("Sample_", 1:n_reps)

dt_out <- data.table::as.data.table(random20_matrix)
data.table::fwrite(dt_out, file = 'luad_random20_markers_100rep.csv')
