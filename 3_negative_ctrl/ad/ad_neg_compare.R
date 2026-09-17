# ! Generate neg ctrl markers

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(file.path(usethis::proj_path(), "3_negative_ctrl/ad"))

set.seed(123)

library(dplyr)
library(data.table)

# ! ov-sc-GSE165897

data_path = '/home/data/sigbridger/benchmark_data/ad'
seurat <- qs::qread(
  file.path(data_path, 'ad_GSE138852_seurat.qs'),
  nthreads = 4
)

# ! bulk
GSE28146 <- qs::qread(
  file.path(data_path, 'ad_bulkdata_GSE28146.qs'),
  nthreads = 4
)
GSE39420 <- qs::qread(
  file.path(data_path, 'ad_bulkdata_GSE39420.qs'),
  nthreads = 4
)
GSE109887 <- qs::qread(
  file.path(data_path, 'ad_bulkdata_GSE109887.qs'),
  nthreads = 4
)


# ! random 20 gene
gene_sc <- rownames(seurat)
gene_GSE28146 <- rownames(GSE28146)
gene_GSE39420 <- rownames(GSE39420)
gene_GSE109887 <- rownames(GSE109887)

common_gene <- Reduce(
  intersect,
  list(gene_sc, gene_GSE28146, gene_GSE39420, gene_GSE109887)
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
data.table::fwrite(dt_out, file = 'ad_random20_markers_100rep.csv')
