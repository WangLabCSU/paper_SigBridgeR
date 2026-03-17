library(dplyr)
library(data.table)


# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/survival/brca/her2/GSE162228"
)

data_path = "/home/data/sigbridger/benchmark_data/brca"
devtools::document("~/R/Project/R_code/SigBridgeR")

# * load data
seurat = qs::qread(file.path(data_path, "seurat_her2.qs"))

bulk = qs::qread(
  file.path(data_path, "brca_bulkdata_GSE162228.qs")
)

pheno = qs::qread(file.path(
  data_path,
  "brca_pheno_GSE162228.qs"
))

surv_data = pheno %>%
  select("overall survival (years):ch1", "characteristics_ch1.5") %>%
  rename("time" := 1) %>%
  mutate(
    status = case_when(
      characteristics_ch1.5 == "alive: Death" ~ 1,
      characteristics_ch1.5 == "alive: Alive" ~ 0
    )
  ) %>%
  select(-"characteristics_ch1.5") %>%
  mutate_all(~ as.numeric(.))

# bulk = bulk[, rownames(surv_data)]

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

seurat_tum = readRDS(
  '/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_HER2Tum.rds'
)
tumor_cells = rownames(seurat_tum[[]])
pos_ratio$benchmark = ifelse(colnames(seurat) %chin% tumor_cells, TRUE, FALSE)
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
# 2349546
# ! ℹ [2025/11/25 08:55:21] Scissor start...
# ! ℹ [2025/11/25 08:55:21] Start from raw data...
# ! ℹ Using "RNA_snn" graph for network.
# ! Warning in asMethod(object) :
# !   sparse->dense coercion: allocating vector of size 7.6 GiB
# ! ℹ [2025/11/25 08:55:29] Normalizing quantiles of data
# ! Warning in asMethod(object) :
# !   sparse->dense coercion: allocating vector of size 2.7 GiB
# ! ℹ [2025/11/25 08:56:26] Subsetting data
# ! ℹ [2025/11/25 08:56:32] Calculating correlation
# ! -------------------------------------------------------------------------------
# ! Five-number summary of correlations:
# ! 0.406719 0.449469 0.461075 0.471802 0.500627
# ! -------------------------------------------------------------------------------
# ! ℹ [2025/11/25 08:57:57] Perform cox regression on the given clinical outcomes...
# ! ✔ [2025/11/25 08:59:32] Statistics data saved to Scissor_inputs.RData.
# ! ℹ [2025/11/25 08:59:33] Screening...

# ! ── At alpha = 0.005 ──

# ! Scissor identified 0 Scissor+ cells and 0 Scissor- cells.
# ! The percentage of selected cell is: 0%
