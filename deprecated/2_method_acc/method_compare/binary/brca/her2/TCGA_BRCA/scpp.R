library(dplyr)
library(data.table)
library(ScPP)
# library(furrr)
library(BiocParallel)


# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/binary/brca/her2/TCGA_BRCA"
)

data_path = "/home/data/sigbridger/benchmark_data/brca"
devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "seurat_her2.qs"))

bulk = qs::qread(
  file.path(data_path, "brca_bulkdata_TCGA.qs")
)

pheno = qs::qread(
  file.path(data_path, "brca_pheno_TCGA.qs"),
  nthreads = 4
)


pheno_bi = mutate(
  pheno,
  sample_type = substr(pheno$sample, 14, 15)
) %>%
  select(sample, sample_type) %>%
  filter(sample_type %in% c("01", "11")) %>%
  mutate(sample_type = ifelse(sample_type == "01", 1, 0))

pheno_bi = setNames(pheno_bi$sample_type, pheno_bi$sample)

cm_samples = intersect(names(pheno_bi), colnames(bulk))
bulk = bulk[, cm_samples]
pheno_bi = pheno_bi[cm_samples]


if (!all(colnames(bulk) == names(pheno_bi))) {
  stop("bulk and pheno_bi not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}

# --------------------------------------------------------------------------------------------
