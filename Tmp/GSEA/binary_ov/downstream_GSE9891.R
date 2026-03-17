# setwd(.rs.api.getActiveDocumentContext()$path |> dirname())
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/GSEA/binary_ov")

data_path = "/home/data/sigbridger/GSEA/ov"

library(irGSEA)

irgsea_score = qs::qread(
  file.path(data_path, "ov_GSE140082_irGSEA_score.qs"),
  nthreads = 8L
)

# ! BULK- GSE9891
# ! SC- GSE165897 - ov
# ! phenotype - binary variable

ov_GSE9891_merged <- qs::qread(
  "/home/data/sigbridger/benchmark_binary/ov/GSE165897/GSE9891_hgsoc_merged_seurat.qs",
  nthreads = 8L
)
meta = ov_GSE9891_merged[[]]

irgsea_score$scissor = meta$scissor
irgsea_score$scPAS = meta$scPAS
irgsea_score$scAB = meta$scAB
irgsea_score$scPAS_RS = meta$scPAS_RS
irgsea_score$scPAS_NRS = meta$scPAS_NRS
irgsea_score$scPAS_Pvalue = meta$scPAS_Pvalue
irgsea_score$scPAS_FDR = meta$scPAS_FDR
irgsea_score$scAB_Subset1 = meta$scAB_Subset1
irgsea_score$scAB_Subset2 = meta$scAB_Subset2
irgsea_score$Subset1_loading = meta$Subset1_loading
irgsea_score$Subset2_loading = meta$Subset2_loading
irgsea_score$scpp = meta$scpp
# irgsea_score$scPP = NULL

labeled_irgsea = SigBridgeR::MergeResult(irgsea_score, ov_GSE9891_merged)

# Wlicox test is perform to all enrichment score matrixes and gene sets
# with adjusted p value &lt; 0.05 are used to integrated through RRA.
# Among them, Gene sets with p value &lt; 0.05 are statistically
# significant and common differential in all gene sets enrichment analysis
# methods. All results are saved in a list.
scissor.dge <- irGSEA::irGSEA.integrate(
  object = labeled_irgsea,
  group.by = "scissor",
  metadata = NULL,
  col.name = NULL,
  method = c(
    "AUCell",
    "UCell",
    "singscore",
    "ssgsea"
    #,         "JASMINE",
    # "viper"
  )
)

cli::cli_alert_success("Done!")

scpas.dge <- irGSEA::irGSEA.integrate(
  object = labeled_irgsea,
  group.by = "scPAS",
  metadata = NULL,
  col.name = NULL,
  method = c(
    "AUCell",
    "UCell",
    "singscore",
    "ssgsea"
    #, "JASMINE",
    # "viper"
  )
)

cli::cli_alert_success("Done!")

scab.dge <- irGSEA::irGSEA.integrate(
  object = labeled_irgsea,
  group.by = "scAB",
  metadata = NULL,
  col.name = NULL,
  method = c(
    "AUCell",
    "UCell",
    "singscore",
    "ssgsea"
    # ,        "JASMINE",
    # "viper"
  )
)

cli::cli_alert_success("Done!")

scpp.dge <- irGSEA::irGSEA.integrate(
  object = labeled_irgsea,
  group.by = "scPP",
  metadata = NULL,
  col.name = NULL,
  method = c(
    "AUCell",
    "UCell",
    "singscore",
    "ssgsea"
    #, "JASMINE",
    # "viper"
  )
)
cli::cli_alert_success("Done!")

output_dir = "/home/data/sigbridger/GSEA/binary_ov"

qs::qsave(
  list(
    scissor = scissor.dge,
    scab = scab.dge,
    scpas = scpas.dge,
    scpp = scpp.dge
  ),
  file = file.path(output_dir, "binary_ov_GSE9891_dge.qs"),
  nthreads = 4L
)
# PID='995694'

cli::cli_h1("All Done!")
