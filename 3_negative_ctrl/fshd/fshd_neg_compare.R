# ! Generate neg ctrl markers

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(file.path(usethis::proj_path(), "3_negative_ctrl/fshd"))

set.seed(123)

library(dplyr)
library(data.table)

# ! ov-sc-GSE122873

data_path = '/home/data/sigbridger/benchmark_data/fshd'
seurat <- qs::qread(
  file.path(data_path, 'fshd_GSE122873_seurat.qs'),
  nthreads = 4
)

# ! bulk
GSE56787 <- qs::qread(
  file.path(data_path, 'fshd_bulkdata_GSE56787.qs'),
  nthreads = 4
)
GSE115650 <- qs::qread(
  file.path(data_path, 'fshd_bulkdata_GSE115650.qs'),
  nthreads = 4
)
GSE140261 <- qs::qread(
  file.path(data_path, 'fshd_bulkdata_GSE140261.qs'),
  nthreads = 4
)


# ! random 20 gene
gene_sc <- rownames(seurat)
gene_GSE56787 <- rownames(GSE56787)
gene_GSE115650 <- rownames(GSE115650)
gene_GSE140261 <- rownames(GSE140261)

common_gene <- Reduce(
  intersect,
  list(gene_sc, gene_GSE56787, gene_GSE115650, gene_GSE140261)
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
data.table::fwrite(dt_out, file = 'fshd_random20_markers_100rep.csv')
