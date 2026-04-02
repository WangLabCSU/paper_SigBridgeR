# ! LUAD bulk-UMI count GSE3141, scRNA-seq data- GSE123902
# ! binary
# ! cannot access the binary phenotype from the article, deprecated

library(dplyr)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_binary/lung/GSE3141")
data_path = "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_data/lung"
devtools::document("~/R/Project/R_code/SigBridgeR")


# Load bulk data & phenotype
bulkdata = qs::qread(file.path(data_path, "lung_bulkdata_GSE3141.qs"))

pheno = data.table::fread(
  "/home/data/data-resource/single-cell/Lung_Cancer/LUAD/lung.cancer.adeno.gse3141.hgu133plus2_entrezcdf.tsv"
)

# check
# bulkdata_ok = BulkPreProcess(
#     bulkdata,
#     gene_symbol_conversion = FALSE,
#     check = FALSE,
#     min_count_threshold = 10,
#     min_gene_expressed = 3,
#     min_total_reads = 1e+05,
#     min_genes_detected = 10000,
#     min_correlation = 0.8,
#     n_top_genes = 500,
#     show_plot_results = TRUE,
#     verbose = TRUE
# )
bulkdata_ok = bulkdata

pheno_ok = pheno %>%
  filter(`Array` %in% colnames(bulkdata_ok)) %>%
  tibble::column_to_rownames("Array") %>%
  select("OS_Time", "OS_Status") %>%
  rename("time" = "OS_Time", "status" = "OS_Status")


# load single-cell data
scdata_ok = qs::qread(
  file.path(data_path, "luad_GSE123902_seurat.qs"),
  nthreads = 4
)

# here `anndata` installed
# reticulate::use_condaenv("base")

# scdata = anndata::read_h5ad(
#     "/home/data/data-resource/single-cell/NC-Zenodo10651059/atlas_dataset/lu_adc_GSE123902.h5ad"
# )
# scdata_ok = SCPreProcess(
#     scdata,
#     quality_control.pattern = '^MT-',
#     dims = 1:40
# )

# qs::qsave(scdata_ok, "../luad_GSE123902_seurat.qs", nthreads = 4)

scissor_result = Screen(
  bulkdata_ok,
  scdata_ok,
  pheno_ok,
  label_type = "scissor_survival",
  phenotype_class = "survival",
  screen_method = "Scissor"
)


qs::qsave(scissor_result, "luad_GSE3141_scissor.qs", nthreads = 4)

scpas_result = Screen(
  bulkdata_ok,
  scdata_ok,
  pheno_ok,
  label_type = "scpas_survival",
  phenotype_class = "survival",
  screen_method = "scPAS"
)

qs::qsave(scpas_result, "luad_GSE3141_scpas.qs", nthreads = 4)

scab_result = Screen(
  bulkdata_ok,
  scdata_ok,
  pheno_ok,
  label_type = "scab_survival",
  phenotype_class = "survival",
  screen_method = "scAB"
)

qs::qsave(scab_result, "luad_GSE3141_scab.qs", nthreads = 4)

scpp_result = Screen(
  bulkdata_ok,
  scdata_ok,
  pheno_ok,
  label_type = "scpp_survival",
  phenotype_class = "survival",
  screen_method = "scPP"
)

qs::qsave(scpp_result, "luad_GSE3141_scpp.qs", nthreads = 4)
