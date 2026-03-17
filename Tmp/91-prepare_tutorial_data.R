setwd("/home/yyx/R/Project/R_code/SigBridgeR/vignettes/example_data")
library(Matrix)
binary_data = readRDS("binary_example_data.rds")

binary_sc = binary_data[[1]]
bianry_bulk = binary_data[[2]]
binary_pheno = binary_data[[3]]

class(binary_sc)
binary_sc_mat = as.matrix(binary_sc)
binary_sc_sparse = Matrix(
  binary_sc_mat
)
class(binary_sc_sparse)


class(bianry_bulk)
bianry_bulk_mat = as.matrix(bianry_bulk)
bianry_bulk_sparse = Matrix(
  bianry_bulk_mat
)

class(bianry_bulk_sparse)

compressed_data = list(binary_sc_mat, bianry_bulk_mat, binary_pheno)

saveRDS(compressed_data, "compressed_binary_example_data.rds")
