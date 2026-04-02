# ! sc - GSE161529
# ! bulk - GSE42568
# ! bulk - GSE162228
# ! binary
# ! TNBC

library(Seurat)
library(dplyr)

devtools::document("/home/yyx/R/Project/R_code/SigBridgeR")

setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_binary/brca/TNBC")

data_path = "/home/data/sigbridger/benchmark_data/brca"
save_path = "/home/data/sigbridger/benchmark_binary/brca/TNBC"

# ? ---- Single cell data preparation ----

# ! BRCA GSE161529-sc
seurat_tnbc = qs::qread(file.path(data_path, "seurat_tnbc.qs"), nthreads = 4)

# # ? ---- Bulk data preparation ----

# * bulk data
brca_bulkdata_GSE42568 = qs::qread(
  file.path(data_path, "brca_bulkdata_GSE42568.qs"),
  nthreads = 4
)

# > dim(brca_bulkdata_GSE42568 )
# [1] 54675   121

# * phenotype
brca_pheno_GSE42568 = qs::qread(
  file.path(data_path, "brca_pheno_GSE42568.qs"),
  nthreads = 4
)
# * binary variables
brca_bi_GSE42568 = setNames(
  case_when(
    brca_pheno_GSE42568$`tissue:ch1` == "breast cancer" ~ 1,
    brca_pheno_GSE42568$`tissue:ch1` == "normal breast" ~ 0
  ),
  brca_pheno_GSE42568$geo_accession
)

# > brca_pheno_GSE42568$`tissue:ch1` |> table()
# breast cancer normal breast
#           104            17
# table(brca_bi_GSE42568)
# brca_bi_GSE42568
#   0   1
#  17 104

filtered_bulk_GSE42568 = brca_bulkdata_GSE42568[, names(brca_bi_GSE42568)]

# > dim(filtered_bulk_GSE42568)
# [1] 54675   104

# scissor_result = Screen(
#     matched_bulk = filtered_bulk_GSE42568,
#     sc_data = seurat_tnbc,
#     phenotype = brca_bi_GSE42568,
#     label_type = "scissor_binary",
#     phenotype_class = "binary",
#     screen_method = "Scissor",
#     path2save_scissor_inputs = "binary_GSE42568_Scissor_inputs.RData",
# )

# qs::qsave(
#     scissor_result,
#     file = "binary_GSE42568_tnbc_scissor_result.qs",
#     nthreads = 4
# )

# scpas_result = Screen(
#     matched_bulk = filtered_bulk_GSE42568,
#     sc_data = seurat_tnbc,
#     phenotype = brca_bi_GSE42568,
#     label_type = "scPAS_binary",
#     phenotype_class = "binary",
#     screen_method = "scPAS"
# )

# qs::qsave(
#     scpas_result,
#     file = "binary_GSE42568_tnbc_scpas_result.qs",
#     nthreads = 4
# )

# scab_result = Screen(
#     matched_bulk = filtered_bulk_GSE42568,
#     sc_data = seurat_tnbc,
#     phenotype = brca_bi_GSE42568,
#     label_type = "scAB_binary",
#     phenotype_class = "binary",
#     screen_method = "scAB"
# )

# qs::qsave(
#     scab_result,
#     file = file.path(save_path, "binary_GSE42568_tnbc_scab_result.qs"),
#     nthreads = 4
# )

# scpp_result = Screen(
#     matched_bulk = filtered_bulk_GSE42568,
#     sc_data = seurat_tnbc,
#     phenotype = brca_bi_GSE42568,
#     label_type = "scPP_binary",
#     phenotype_class = "binary",
#     screen_method = "scPP"
# )

# qs::qsave(
#     scpp_result,
#     file = "binary_GSE42568_tnbc_scpp_result.qs",
#     nthreads = 4
# )

for (method in c(
  "Scissor",
  "scAB",
  "scPAS",
  "scPP",
  "DEGAS",
  "PIPET",
  "LP_SGL"
)) {
  rlang::try_fetch(
    {
      res <- Screen(
        filtered_bulk_GSE42568,
        seurat_tnbc,
        brca_bi_GSE42568,
        label_type = paste0(method, "_binary"),
        phenotype_class = "binary",
        screen_method = method,
        path2save_scissor_inputs = NULL
      )

      qs::qsave(
        res,
        file.path(
          data_path,
          "TNBC",
          paste0("binary_GSE42568_tnbc_", tolower(method), "_result.qs")
        ),
        nthreads = 4
      )
    },
    error = function(e) {
      print(e$message)
      cli::cli_h1("{method} failed")
    }
  )
}

# ? ---- Another bulk data ----

brca_bulkdata_GSE162228 = qs::qread(
  file.path(data_path, "brca_bulkdata_GSE162228.qs"),
  nthreads = 4
)

brca_pheno_GSE162228 = qs::qread(
  file.path(data_path, "brca_pheno_GSE162228.qs"),
  nthreads = 4
)
brca_bi_GSE162228 = setNames(
  case_when(
    brca_pheno_GSE162228$`relapse status:ch1` == "relapse" ~ 1,
    brca_pheno_GSE162228$`relapse status:ch1` == "non-relapse" ~ 0
  ),
  brca_pheno_GSE162228$geo_accession
)
# brca_bi_GSE162228 |> table()
# # brca_bi_GSE162228
# #   0   1
# # 111  22

# scissor_result = Screen(
#     matched_bulk = brca_bulkdata_GSE162228,
#     sc_data = seurat_tnbc,
#     phenotype = brca_bi_GSE162228,
#     label_type = "scissor_binary",
#     phenotype_class = "binary",
#     screen_method = "Scissor",
#     path2save_scissor_inputs = "binary_GSE162228_Scissor_inputs.RData",
# )

# qs::qsave(
#     scissor_result,
#     file = "binary_GSE162228_tnbc_scissor_result.qs",
#     nthreads = 4
# )

# scpas_result = Screen(
#     matched_bulk = brca_bulkdata_GSE162228,
#     sc_data = seurat_tnbc,
#     phenotype = brca_bi_GSE162228,
#     label_type = "scPAS_binary",
#     phenotype_class = "binary",
#     screen_method = "scPAS"
# )

# qs::qsave(
#     scpas_result,
#     file = "binary_GSE162228_tnbc_scpas_result.qs",
#     nthreads = 4
# )

# scab_result = Screen(
#     matched_bulk = brca_bulkdata_GSE162228,
#     sc_data = seurat_tnbc,
#     phenotype = brca_bi_GSE162228,
#     label_type = "scAB_binary",
#     phenotype_class = "binary",
#     screen_method = "scAB"
# )

# qs::qsave(
#     scab_result,
#     file = file.path(save_path, "binary_GSE162228_tnbc_scab_result.qs"),
#     nthreads = 4
# )

# scpp_result = Screen(
#     matched_bulk = brca_bulkdata_GSE162228,
#     sc_data = seurat_tnbc,
#     phenotype = brca_bi_GSE162228,
#     label_type = "scPP_binary",
#     phenotype_class = "binary",
#     screen_method = "scPP"
# )

# qs::qsave(
#     scpp_result,
#     file = "binray_GSE162228_tnbc_scpp_result.qs",
#     nthreads = 4
# )

# PID=2543637

for (method in c(
  "Scissor",
  "scAB",
  "scPAS",
  "scPP",
  "DEGAS",
  "PIPET",
  "LP_SGL"
)) {
  rlang::try_fetch(
    {
      res <- Screen(
        brca_bulkdata_GSE162228,
        seurat_tnbc,
        brca_bi_GSE162228,
        label_type = paste0(method, "_binary"),
        phenotype_class = "binary",
        screen_method = method,
        path2save_scissor_inputs = NULL
      )

      qs::qsave(
        res,
        file.path(
          data_path,
          "TNBC",
          paste0("binary_GSE162228_tnbc_", tolower(method), "_result.qs")
        ),
        nthreads = 4
      )
    },
    error = function(e) {
      print(e$message)
      cli::cli_h1("{method} failed")
    }
  )
}


# ? ---- Another bulk data ----

tcga_bulkdata = qs::qread(
  file.path(data_path, "brca_bulkdata_TCGA.qs"),
  nthreads = 4
)
tcga_pheno = qs::qread(
  file.path(data_path, "brca_pheno_TCGA.qs"),
  nthreads = 4
)
cm_samples = intersect(tcga_pheno$sample, colnames(tcga_bulkdata))

brca_bi_tcga = mutate(
  tcga_pheno,
  sample_type = substr(tcga_pheno$sample, 14, 15)
) %>%
  select(sample, sample_type) %>%
  filter(sample_type %in% c("01", "11"), sample %in% cm_samples) %>%
  mutate(sample_type = ifelse(sample_type == "01", 1, 0))
brca_bi_tcga = setNames(brca_bi_tcga$sample_type, brca_bi_tcga$sample)

tcga_bulkdata = tcga_bulkdata[, names(brca_bi_tcga)]

# table(brca_bi_tcga)
# # brca_bi_tcga
# #    0    1
# #  113 1089

# scissor_result = Screen(
#     matched_bulk = tcga_bulkdata,
#     sc_data = seurat_tnbc,
#     phenotype = brca_bi_tcga,
#     label_type = "binary_scissor",
#     phenotype_class = "binary",
#     screen_method = "Scissor",
#     path2save_scissor_inputs = "binary_tcga_Scissor_inputs.RData",
# )

# qs::qsave(
#     scissor_result,
#     file = file.path(save_path, "binary_tnbc_tcga_scissor_result.qs"),
#     nthreads = 4
# )

# scpas_result = Screen(
#     matched_bulk = tcga_bulkdata,
#     sc_data = seurat_tnbc,
#     phenotype = brca_bi_tcga,
#     label_type = "scPAS_binary",
#     phenotype_class = "binary",
#     screen_method = "scPAS"
# )

# qs::qsave(
#     scpas_result,
#     file = file.path(save_path, "binary_tnbc_tcga_scpas_result.qs"),
#     nthreads = 4
# )

# scab_result = Screen(
#     matched_bulk = tcga_bulkdata,
#     sc_data = seurat_tnbc,
#     phenotype = brca_bi_tcga,
#     label_type = "scAB_binary",
#     phenotype_class = "binary",
#     screen_method = "scAB"
# )

# qs::qsave(
#     scab_result,
#     file = file.path(save_path, "binary_tnbc_tcga_scab_result.qs"),
#     nthreads = 4
# )

# scpp_result = Screen(
#   matched_bulk = tcga_bulkdata,
#   sc_data = seurat_tnbc,
#   phenotype = brca_bi_tcga,
#   label_type = "scPP_binary",
#   phenotype_class = "binary",
#   screen_method = "scPP"
# )

# qs::qsave(
#   scpp_result,
#   file = file.path(save_path, "binary_tnbc_tcga_scpp_result.qs"),
#   nthreads = 4
# )
# # 266009

for (method in c(
  "Scissor",
  "scAB",
  "scPAS",
  "scPP",
  "DEGAS",
  "PIPET",
  "LP_SGL"
)) {
  rlang::try_fetch(
    {
      res <- Screen(
        tcga_bulkdata,
        seurat_tnbc,
        brca_bi_tcga,
        label_type = paste0(method, "_binary"),
        phenotype_class = "binary",
        screen_method = method,
        path2save_scissor_inputs = NULL
      )

      qs::qsave(
        res,
        file.path(
          data_path,
          "TNBC",
          paste0("binary_tcga_tnbc_", tolower(method), "_result.qs")
        ),
        nthreads = 4
      )
    },
    error = function(e) {
      print(e$message)
      cli::cli_h1("{method} failed")
    }
  )
}
