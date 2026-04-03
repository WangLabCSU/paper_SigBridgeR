# ! GSE165897-sc
# ! GSE9891-GSE140082-bulk, GSE32062没有合适的二元变量数据

library(Seurat)
library(dplyr)

library(SigBridgeR)

setwd("~/R/Project/R_code/SigBridgeR/Tmp/benchmark_binary/ov/GSE165897")

data_path = "/home/data/sigbridger/benchmark_data/ov"
save_path = "/home/data/sigbridger/benchmark_binary/ov/GSE165897"
# *sc
GSE165897_seurat = qs::qread(
  file.path(data_path, "hgsoc_GSE165897_seurat.qs"),
  nthreads = 4
)


# *pheno
GSE9891_pheno_raw = qs::qread(
  file.path(data_path, "ov_pheno_GSE9891.qs"),
  nthreads = 4
)

# *binary data for screen
GSE9891_pheno = setNames(
  case_when(
    GSE9891_pheno_raw$characteristics_ch1.1 == "Type : LMP" ~ 0,
    GSE9891_pheno_raw$characteristics_ch1.1 == "Type : Malignant" ~ 1
  ),
  GSE9891_pheno_raw$geo_accession
)
# > table(GSE9891_pheno_raw$characteristics_ch1.1)
#       Type : LMP Type : Malignant
#               18              267
# > table(GSE9891_pheno)
# GSE9891_pheno
#   0   1
#  18 267

# *bulk
GSE9891_bulk = qs::qread(
  file.path(data_path, "ov_bulkdata_GSE9891.qs"),
  nthreads = 4
)[,
  names(GSE9891_pheno)
]

# * bulk
# ! GSE32062表型没有好的二元变量数据
# GSE32062_bulk = qs::qread(
#     "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_data/ov/ov_bulkdata_GSE32062_GPL6480.qs",
#     nthreads = 4
# )
GSE140082_bulk = qs::qread(
  file.path(data_path, "ov_bulkdata_GSE140082.qs"),
  nthreads = 4
)
# GSM4153781 has NA values, find median and replace NA

# # * pheno
# # ! GSE32062表型没有好的二元变量数据
# GSE32062_pheno = qs::qread(
#     "/home/data/sigbridger/benchmark_data/ov/ov_pheno_GSE32062.qs",
#     nthreads = 4
# )

GSE140082_pheno_raw = qs::qread(
  file.path(data_path, "ov_pheno_GSE140082.qs"),
  nthreads = 4
)
GSE140082_pheno_raw = GSE140082_pheno_raw[
  GSE140082_pheno_raw$`newgrade:ch1` != "NA" &
    !is.na(GSE140082_pheno_raw$`newgrade:ch1`),
]
GSE140082_pheno = setNames(
  case_when(
    GSE140082_pheno_raw$`newgrade:ch1` == "high.grade" ~ 1,
    GSE140082_pheno_raw$`newgrade:ch1` == "low.grade" ~ 0
  ),
  GSE140082_pheno_raw$geo_accession
)

# ! 匹配
GSE140082_bulk = GSE140082_bulk[, names(GSE140082_pheno)]

# GSM4153781 has NA values, find median and replace NA
median = median(GSE140082_bulk[, "GSM4153781"], na.rm = TRUE)
na_indices <- which(is.na(GSE140082_bulk[, "GSM4153781"]))
GSE140082_bulk[na_indices, "GSM4153781"] <- median

# --------------------------------------------------------------------

# # * screen
# scissor_result = SigBridgeR::Screen(
#     matched_bulk = GSE9891_bulk,
#     sc_data = GSE165897_seurat,
#     phenotype = GSE9891_pheno,
#     label_type = "Malignant_or_LMP",
#     phenotype_class = "binary",
#     screen_method = "Scissor",
#     path2save_scissor_inputs = "GSE165897_GSE9891_Scissor_inputs.RData",
# )

# qs::qsave(scissor_result, "GSE165897_GSE9891_Scissor_result.qs")

# scissor_result = SigBridgeR::Screen(
#     matched_bulk = GSE32062_bulk,
#     sc_data = GSE165897_seurat,
#     phenotype = GSE32062_pheno,
#     label_type = "scissor",
#     phenotype_class = "binary",
#     screen_method = "Scissor",
#     path2save_scissor_inputs = "GSE165897_GSE32062_Scissor_inputs.RData",
# )

# qs::qsave(scissor_result, "GSE165897_GSE32062_Scissor_result.qs")

# scissor_result = SigBridgeR::Screen(
#     matched_bulk = GSE140082_bulk,
#     sc_data = GSE165897_seurat,
#     phenotype = GSE140082_pheno,
#     label_type = "bevacizumab_or_standard",
#     phenotype_class = "binary",
#     screen_method = "Scissor",
#     path2save_scissor_inputs = "GSE165897_GSE140082_Scissor_inputs.RData",
# )

# qs::qsave(scissor_result, "GSE165897_GSE140082_Scissor_result.qs")

# # # *scPAS
# scpas_result = Screen(
#     matched_bulk = GSE9891_bulk,
#     sc_data = GSE165897_seurat,
#     phenotype = GSE9891_pheno,
#     label_type = "Malignant_or_LMP",
#     phenotype_class = "binary",
#     screen_method = "scPAS"
# )

# qs::qsave(
#     scpas_result,
#     file = "GSE165897_GSE9891_scPAS_result.qs",
#     nthreads = 4
# )

# scpas_result = Screen(
#     matched_bulk = GSE32062_bulk,
#     sc_data = GSE165897_seurat,
#     phenotype = GSE32062_pheno,
#     label_type = "scPAS",
#     phenotype_class = "binary",
#     screen_method = "scPAS"
# )

# qs::qsave(
#     scpas_result,
#     file = "GSE165897_GSE32062_scPAS_result.qs",
#     nthreads = 4
# )

# scpas_result = Screen(
#     matched_bulk = GSE140082_bulk,
#     sc_data = GSE165897_seurat,
#     phenotype = GSE140082_pheno,
#     label_type = "bevacizumab_or_standard",
#     phenotype_class = "binary",
#     screen_method = "scPAS"
# )

# qs::qsave(
#     scpas_result,
#     file = "GSE165897_GSE140082_scPAS_result.qs",
#     nthreads = 4
# )

# *scAB
# scab_result = Screen(
#   matched_bulk = GSE9891_bulk,
#   sc_data = GSE165897_seurat,
#   phenotype = GSE9891_pheno,
#   label_type = "Malignant_or_LMP",
#   phenotype_class = "binary",
#   screen_method = "scAB"
# )

# qs::qsave(
#   scab_result,
#   file = file.path(save_path, "GSE165897_GSE9891_scAB_result.qs"),
#   nthreads = 4
# )

# scab_result = Screen(
#   matched_bulk = GSE140082_bulk,
#   sc_data = GSE165897_seurat,
#   phenotype = GSE140082_pheno,
#   label_type = "bevacizumab_or_standard",
#   phenotype_class = "binary",
#   screen_method = "scAB"
# )

# qs::qsave(
#   scab_result,
#   file = file.path(save_path, "GSE165897_GSE140082_scAB_result.qs"),
#   nthreads = 4
# )

# # *scPP
# scpp_result = Screen(
#     matched_bulk = GSE9891_bulk,
#     sc_data = GSE165897_seurat,
#     phenotype = GSE9891_pheno,
#     label_type = "Malignant_or_LMP",
#     phenotype_class = "binary",
#     screen_method = "scPP"
# )

# qs::qsave(
#     scpp_result,
#     file = "GSE165897_GSE9891_scPP_result.qs",
#     nthreads = 4
# )

# scpp_result = Screen(
#     matched_bulk = GSE140082_bulk,
#     sc_data = GSE165897_seurat,
#     phenotype = GSE140082_pheno,
#     label_type = "bevacizumab_or_standard",
#     phenotype_class = "binary",
#     screen_method = "scPP"
# )

# qs::qsave(
#     scpp_result,
#     file = "GSE165897_GSE140082_scpp_result.qs",
#     nthreads = 4
# )

# # ! scPP is not supported in GSE32062
# scpp_result = Screen(
#     matched_bulk = GSE32062_bulk,
#     sc_data = GSE165897_seurat,
#     phenotype = GSE32062_pheno,
#     label_type = "scPP",
#     phenotype_class = "binary",
#     screen_method = "scPP"
# )

# qs::qsave(
#     scpp_result,
#     file = "GSE165897_GSE32062_scPP_result.qs",
#     nthreads = 4
# )

# PID=283845

SigBridgeR::setThreads(
  4L,
  tf_config = list(xla = TRUE, inter_op = 4L, intra_op = 4L)
)

for (method in c(
  "Scissor",
  "scAB",
  #   "scPAS",
  "scPP",
  "DEGAS",
  #   "PIPET",
  "LP_SGL"
)) {
  #   rlang::try_fetch(
  #     {
  #       screen_9891 = Screen(
  #         matched_bulk = GSE9891_bulk,
  #         sc_data = GSE165897_seurat,
  #         phenotype = GSE9891_pheno,
  #         label_type = method,
  #         phenotype_class = "binary",
  #         screen_method = method
  #       )
  #       qs::qsave(
  #         screen_9891,
  #         file = file.path(
  #           save_path,
  #           paste0("GSE165897_GSE9891_", method, "_result.qs")
  #         ),
  #         nthreads = 4
  #       )
  #     },
  #     error = function(e) {
  #       cli::cli_h1("{method} failed")
  #     }
  #   )

  rlang::try_fetch(
    {
      screen_140082 = Screen(
        matched_bulk = GSE140082_bulk,
        sc_data = GSE165897_seurat,
        phenotype = GSE140082_pheno,
        label_type = method,
        phenotype_class = "binary",
        screen_method = method
      )

      qs::qsave(
        screen_140082,
        file = file.path(
          save_path,
          paste0("GSE165897_GSE140082_", method, "_result.qs")
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
# PID = 2588004
