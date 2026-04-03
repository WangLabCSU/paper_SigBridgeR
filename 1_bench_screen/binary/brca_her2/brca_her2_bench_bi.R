# ! sc - GSE161529
# ! bulk - GSE42568 , GSE162228 , TCGA_BRCA
# ! binary
# ! 20250910

library(Seurat)
library(dplyr)
library(SigBridgeR)

setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_binary/brca/HER2")

data_path <- "/home/data/sigbridger/benchmark_data/brca"
save_path <- "/home/data/sigbridger/benchmark_binary/brca/HER2"

# ? ---- Single cell data preparation ----

# ! BRCA HER2 GSE161529-sc
seurat_her2 <- qs::qread(file.path(data_path, "seurat_her2.qs"), nthreads = 4)

# ? ---- Bulk data preparation ----

# *bulk data
brca_bulkdata_GSE42568 <- qs::qread(
  file.path(data_path, "brca_bulkdata_GSE42568.qs"),
  nthreads = 4
)

# > dim(brca_bulkdata_GSE42568 )
# [1] 54675   121

# *phenotype
brca_pheno_GSE42568 <- qs::qread(
  file.path(data_path, "brca_pheno_GSE42568.qs"),
  nthreads = 4
)

# * binary variables
brca_bi_GSE42568 <- setNames(
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

filtered_bulk_GSE42568 <- brca_bulkdata_GSE42568[, names(brca_bi_GSE42568)]

# > dim(filtered_bulk_GSE42568)
# [1] 54675   104

# scissor_result = Screen(
#     matched_bulk = filtered_bulk_GSE42568,
#     sc_data = seurat_her2,
#     phenotype = brca_bi_GSE42568,
#     label_type = "binary_scissor",
#     phenotype_class = "binary",
#     screen_method = "Scissor",
#     path2save_scissor_inputs = "binary_GSE42569_Scissor_inputs.RData",
# )

# qs::qsave(
#     scissor_result,
#     file = "binary_her2_GSE42568_scissor_result.qs",
#     nthreads = 4
# )

# scpas_result = Screen(
#     matched_bulk = filtered_bulk_GSE42568,
#     sc_data = seurat_her2,
#     phenotype = brca_bi_GSE42568,
#     label_type = "binary_scPAS",
#     phenotype_class = "binary",
#     screen_method = "scPAS"
# )

# qs::qsave(
#     scpas_result,
#     file = "binary_her2_GSE42568_scpas_result.qs",
#     nthreads = 4
# )

# scab_result = Screen(
#     matched_bulk = filtered_bulk_GSE42568,
#     sc_data = seurat_her2,
#     phenotype = brca_bi_GSE42568,
#     label_type = "binary_scAB",
#     phenotype_class = "binary",
#     screen_method = "scAB"
# )

# qs::qsave(
#     scab_result,
#     file = file.path(save_path, "binary_her2_GSE42568_scab_result.qs"),
#     nthreads = 4
# )

# scpp_result = Screen(
#     matched_bulk = filtered_bulk_GSE42568,
#     sc_data = seurat_her2,
#     phenotype = brca_bi_GSE42568,
#     label_type = "binary_scPP",
#     phenotype_class = "binary",
#     screen_method = "scPP"
# )

# qs::qsave(
#     scpp_result,
#     file = "binary_her2_GSE42568_scpp_result.qs",
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
        seurat_her2,
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
          paste0("binary_her2_GSE42568_", tolower(method), "_result.qs")
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

brca_bulkdata_GSE162228 <- qs::qread(
  file.path(data_path, "brca_bulkdata_GSE162228.qs"),
  nthreads = 4
)

brca_pheno_GSE162228 <- qs::qread(
  file.path(data_path, "brca_pheno_GSE162228.qs"),
  nthreads = 4
)
brca_bi_GSE162228 <- setNames(
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
#     sc_data = seurat_her2,
#     phenotype = brca_bi_GSE162228,
#     label_type = "binary_scissor",
#     phenotype_class = "binary",
#     screen_method = "Scissor",
#     path2save_scissor_inputs = "binary_GSE162228_Scissor_inputs.RData",
# )

# qs::qsave(
#     scissor_result,
#     file = "binary_her2_GSE162228_scissor_result.qs",
#     nthreads = 4
# )

# scpas_result = Screen(
#     matched_bulk = brca_bulkdata_GSE162228,
#     sc_data = seurat_her2,
#     phenotype = brca_bi_GSE162228,
#     label_type = "scPAS_binary",
#     phenotype_class = "binary",
#     screen_method = "scPAS"
# )

# qs::qsave(
#     scpas_result,
#     file = "binary_her2_GSE162228_scpas_result.qs",
#     nthreads = 4
# )

# scab_result = Screen(
#     matched_bulk = brca_bulkdata_GSE162228,
#     sc_data = seurat_her2,
#     phenotype = brca_bi_GSE162228,
#     label_type = "scAB_binary",
#     phenotype_class = "binary",
#     screen_method = "scAB"
# )

# qs::qsave(
#     scab_result,
#     file = file.path(save_path, "binary_her2_GSE162228_scab_result.qs"),
#     nthreads = 4
# )

# scpp_result = Screen(
#     matched_bulk = brca_bulkdata_GSE162228,
#     sc_data = seurat_her2,
#     phenotype = brca_bi_GSE162228,
#     label_type = "scPP_binary",
#     phenotype_class = "binary",
#     screen_method = "scPP"
# )

# qs::qsave(
#     scpp_result,
#     file = "binary_her2_GSE162228_scpp_result.qs",
#     nthreads = 4
# )

# PID=2533433

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
        seurat_her2,
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
          paste0("binary_her2_GSE162228_", tolower(method), "_result.qs")
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

tcga_bulkdata <- qs::qread(
  file.path(data_path, "brca_bulkdata_TCGA.qs"),
  nthreads = 4
)
tcga_pheno <- qs::qread(
  file.path(data_path, "brca_pheno_TCGA.qs"),
  nthreads = 4
)
cm_samples <- intersect(tcga_pheno$sample, colnames(tcga_bulkdata))

brca_bi_tcga <- mutate(
  tcga_pheno,
  sample_type = substr(tcga_pheno$sample, 14, 15)
) %>%
  select(sample, sample_type) %>%
  filter(sample_type %in% c("01", "11"), sample %in% cm_samples) %>%
  mutate(sample_type = as.integer(sample_type == "01"))
brca_bi_tcga <- setNames(brca_bi_tcga$sample_type, brca_bi_tcga$sample)

tcga_bulkdata <- tcga_bulkdata[, names(brca_bi_tcga)]

# table(brca_bi_tcga)
# # brca_bi_tcga
# #    0    1
# #  113 1089

# scissor_result = Screen(
#     matched_bulk = tcga_bulkdata,
#     sc_data = seurat_her2,
#     phenotype = brca_bi_tcga,
#     label_type = "binary_scissor",
#     phenotype_class = "binary",
#     screen_method = "Scissor",
#     path2save_scissor_inputs = "binary_tcga_Scissor_inputs.RData",
# )

# qs::qsave(
#     scissor_result,
#     file = file.path(save_path, "binary_her2_tcga_scissor_result.qs"),
#     nthreads = 4
# )

# scpas_result = Screen(
#     matched_bulk = tcga_bulkdata,
#     sc_data = seurat_her2,
#     phenotype = brca_bi_tcga,
#     label_type = "scPAS_binary",
#     phenotype_class = "binary",
#     screen_method = "scPAS"
# )

# qs::qsave(
#     scpas_result,
#     file = file.path(save_path, "binary_her2_tcga_scpas_result.qs"),
#     nthreads = 4
# )

# scab_result = Screen(
#     matched_bulk = tcga_bulkdata,
#     sc_data = seurat_her2,
#     phenotype = brca_bi_tcga,
#     label_type = "scAB_binary",
#     phenotype_class = "binary",
#     screen_method = "scAB"
# )

# qs::qsave(
#     scab_result,
#     file = file.path(save_path, "binary_her2_tcga_scab_result.qs"),
#     nthreads = 4
# )

# scpp_result = Screen(
#   matched_bulk = tcga_bulkdata,
#   sc_data = seurat_her2,
#   phenotype = brca_bi_tcga,
#   label_type = "scPP_binary",
#   phenotype_class = "binary",
#   screen_method = "scPP"
# )

# qs::qsave(
#   scpp_result,
#   file = file.path(save_path, "binary_her2_tcga_scpp_result.qs"),
#   nthreads = 4
# )

# cli::cli_alert_success("Benchmark test finished.")
# # 193675
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
        seurat_her2,
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
          paste0("binary_her2_tcga_", tolower(method), "_result.qs")
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
