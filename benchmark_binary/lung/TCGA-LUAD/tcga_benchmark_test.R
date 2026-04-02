# ! LUAD bulk-UMI count TCGA-LUAD, scRNA-seq data- GSE123902
library(dplyr)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_binary/lung/TCGA-LUAD"
)
data_path = "/home/data/sigbridger/benchmark_data/lung"
save_path = "/home/data/sigbridger/benchmark_binary/lung/TCGA-LUAD"
devtools::document("~/R/Project/R_code/SigBridgeR")


# Load bulk data & phenotype
bulkdata_ok = qs::qread(file.path(data_path, "TCGA_LUAD_bulkdata.qs"))

pheno = qs::qread(file.path(data_path, "TCGA_LUAD_pheno.qs"))

# table(sample_type)
# sample_type
#  01  02  11
# 504   2  59

pheno_ok = mutate(pheno, sample_type = substr(pheno$sample, 14, 15)) %>%
  select(sample, sample_type) %>%
  filter(sample_type %in% c("01", "11")) %>%
  mutate(sample_type = ifelse(sample_type == "01", 1, 0))
pheno_ok = setNames(pheno_ok$sample_type, pheno_ok$sample)
# TCGA-05-4244-01 TCGA-05-4250-01 TCGA-05-4382-01 TCGA-05-4384-01 TCGA-05-4389-01 TCGA-05-4390-01
#               1               1               1               1               1               1

bulkdata_ok = bulkdata_ok[, names(pheno_ok)]


# #  check
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

all(colnames(bulkdata_ok) == names(pheno_ok))


# load single-cell data
scdata_ok = qs::qread(
  file.path(data_path, "luad_GSE123902_seurat.qs"),
  nthreads = 4
)

# scissor_result = Screen(
#     bulkdata_ok,
#     scdata_ok,
#     pheno_ok,
#     label_type = "scissor_binary",
#     phenotype_class = "binary",
#     screen_method = "Scissor"
# )

# qs::qsave(scissor_result, "binary_luad_tcga_scissor.qs", nthreads = 4)

# scpas_result = Screen(
#     bulkdata_ok,
#     scdata_ok,
#     pheno_ok,
#     label_type = "scpas_binary",
#     phenotype_class = "binary",
#     screen_method = "scPAS"
# )

# qs::qsave(scpas_result, "binary_luad_tcga_scpas.qs", nthreads = 4)

# scab_result = Screen(
#     bulkdata_ok,
#     scdata_ok,
#     pheno_ok,
#     label_type = "scab_binary",
#     phenotype_class = "binary",
#     screen_method = "scAB"
# )

# qs::qsave(
#     scab_result,
#     file.path(save_path, "binary_luad_tcga_scab.qs"),
#     nthreads = 4
# )

# scpp_result = Screen(
#     bulkdata_ok,
#     scdata_ok,
#     pheno_ok,
#     label_type = "scpp_binary",
#     phenotype_class = "binary",
#     screen_method = "scPP"
# )

# qs::qsave(scpp_result, "binary_luad_tcga_scpp.qs", nthreads = 4)

# PID = 205041

SigBridgeR::setThreads(
  8L,
  tf_config = list(xla = TRUE, intra_op = 8L, inter_op = 8L)
)

for (method in c(
  "Scissor",
  "scAB",
  "scPAS",
  "scPP",
  "DEGAS",
  "LP_SGL",
  "PIPET"
)) {
  rlang::try_fetch(
    {
      res <- Screen(
        bulkdata_ok,
        scdata_ok,
        pheno_ok,
        label_type = paste0(method, "_survival"),
        phenotype_class = "binary",
        screen_method = method,
        path2save_scissor_inputs = NULL
      )

      qs::qsave(
        res,
        file.path(save_path, paste0("luad_TCGA_LUAD_", tolower(method), ".qs")),
        nthreads = 4
      )
    },
    error = function(e) {
      print(e$message)
      cli::cli_h1("{method} failed")
    }
  )
}
