library(dplyr)
library(data.table)


# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/binary/lung/TCGA_LUAD"
)

data_path = "/home/data/sigbridger/benchmark_data/lung"
devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "luad_GSE123902_seurat.qs"))

bulk = qs::qread(
  file.path(data_path, "TCGA_LUAD_bulkdata.qs")
)

pheno = qs::qread(file.path(data_path, "TCGA_LUAD_pheno.qs"))

pheno_bi <- mutate(pheno, sample_type = substr(pheno$sample, 14, 15)) %>%
  select(sample, sample_type) %>%
  filter(sample_type %in% c("01", "11")) %>%
  mutate(sample_type = ifelse(sample_type == "01", 1, 0))
pheno_bi <- setNames(pheno_bi$sample_type, pheno_bi$sample)

bulk = bulk[, names(pheno_bi)]

if (!all(colnames(bulk) == names(pheno_bi))) {
  stop("bulk and pheno_bi not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}
if (anyNA(pheno_bi)) {
  stop("pheno_bi has NA")
}


# ! -------------------- scissor ---------------------

# * random search, 100 times
set.seed(123)


# * run scissor

scissor_result = Screen(
  matched_bulk = bulk,
  sc_data = seurat,
  phenotype = pheno_bi,
  label_type = "tumor",
  phenotype_class = "binary",
  screen_method = "Scissor",
  alpha = NULL,
  path2save_scissor_cache = NULL
)


pos_ratio = (scissor_result$scRNA_data$scissor == "Positive")
pos_ratio = as.data.frame(pos_ratio)
rownames(pos_ratio) = colnames(seurat)

pos_ratio$benchmark = ifelse(seurat$cnv_status == 'tumor', TRUE, FALSE)
source(
  '/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/compute_metrics.R'
)
res = compute_metrics(pos_ratio)

data.table::fwrite(
  res,
  file = "scissor_acc.csv",
  row.names = TRUE,
  col.names = TRUE
)
# 608572
