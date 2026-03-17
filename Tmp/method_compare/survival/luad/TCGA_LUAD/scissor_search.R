library(dplyr)
library(data.table)


# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/survival/luad/TCGA_LUAD"
)

data_path = "/home/data/sigbridger/benchmark_data/lung"
devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "luad_GSE123902_seurat.qs"))

bulk = qs::qread(
  file.path(data_path, "TCGA_LUAD_bulkdata.qs")
)

pheno = qs::qread(file.path(data_path, "TCGA_LUAD_pheno.qs"))

surv_data = select(
  pheno,
  "OS.time",
  "OS"
) %>%
  rename("time" = 1, "status" = 2) %>%
  mutate_all(~ as.numeric(.))
rownames(surv_data) = pheno$sample


if (!all(colnames(bulk) == rownames(surv_data))) {
  stop("bulk and surv_data not match")
}
if (anyNA(bulk)) {
  stop("bulk has NA")
}

# ! -------------------- scissor ---------------------

# * random search, 100 times
set.seed(123)


# * run scissor

scissor_result = Screen(
  path2load_scissor_cache = 'Scissor_inputs.RData',
  matched_bulk = bulk,
  sc_data = seurat,
  phenotype = surv_data,
  label_type = "OS (M)_survival",
  phenotype_class = "survival",
  screen_method = "Scissor",
  alpha = NULL
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
# 3112049
