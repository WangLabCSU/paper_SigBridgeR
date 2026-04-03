# ! TCGA_LUAD

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  file.path(usethis::proj_path(), "2_method_acc/lung")
)
data_path <- "/home/data/sigbridger/benchmark_data/lung"

# * load data
seurat <- qs::qread(file.path(data_path, "luad_GSE123902_seurat.qs"))

bulk <- qs::qread(
  file.path(data_path, "TCGA_LUAD_bulkdata.qs")
)
bulk <- log2(bulk + 1)

pheno <- qs::qread(file.path(data_path, "TCGA_LUAD_pheno.qs"))

pheno_bi <- mutate(pheno, sample_type = substr(pheno$sample, 14, 15)) %>%
  select(sample, sample_type) %>%
  filter(sample_type %in% c("01", "11")) %>%
  mutate(sample_type = as.integer(sample_type == "01"))
pheno_bi <- setNames(pheno_bi$sample_type, pheno_bi$sample)

bulk <- bulk[, names(pheno_bi)]

if (!all(colnames(bulk) == names(pheno_bi))) {
  stop("bulk and pheno_bi not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}
if (anyNA(pheno_bi)) {
  stop("pheno_bi has NA")
}
