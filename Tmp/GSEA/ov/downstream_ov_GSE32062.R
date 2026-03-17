# setwd(.rs.api.getActiveDocumentContext()$path |> dirname())
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/GSEA/ov")

irgsea_score = qs::qread(
  "/home/data/sigbridger/GSEA/ov/ov_GSE140082_irGSEA_score.qs",
  nthreads = 8L
)

ov_GSE32062 = qs::qread(
  "/home/data/sigbridger/benchmark_data/ov/GSE165897/GSE32062_soc_merged_seurat.qs",
  nthreads = 4L
)
meta = ov_GSE32062[[]]

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
# irgsea_score$scpp = meta$scpp # ! 没有scPP
irgsea_score$scPP = NULL

# Wlicox test is perform to all enrichment score matrixes and gene sets
# with adjusted p value &lt; 0.05 are used to integrated through RRA.
# Among them, Gene sets with p value &lt; 0.05 are statistically
# significant and common differential in all gene sets enrichment analysis
# methods. All results are saved in a list.
scissor.dge <- irGSEA::irGSEA.integrate(
  object = irgsea_score,
  group.by = "scissor",
  metadata = NULL,
  col.name = NULL,
  method = c("AUCell", "UCell", "singscore", "ssgsea")
)

cli::cli_alert_success("Done!")

scpas.dge <- irGSEA::irGSEA.integrate(
  object = irgsea_score,
  group.by = "scPAS",
  metadata = NULL,
  col.name = NULL,
  method = c("AUCell", "UCell", "singscore", "ssgsea")
)

cli::cli_alert_success("Done!")

scab.dge <- irGSEA::irGSEA.integrate(
  object = irgsea_score,
  group.by = "scAB",
  metadata = NULL,
  col.name = NULL,
  method = c("AUCell", "UCell", "singscore", "ssgsea")
)

cli::cli_alert_success("Done!")

# scpp.dge <- irGSEA::irGSEA.integrate(
#     object = irgsea_score,
#     group.by = "scPP",
#     metadata = NULL,
#     col.name = NULL,
#     method = c("AUCell", "UCell", "singscore", "ssgsea")
# )
cli::cli_alert_success("Skip scPP due to no Positive cells available!")

qs::qsave(
  list(
    scissor = scissor.dge,
    scab = scab.dge,
    scpas = scpas.dge
    # ,        scpp = scpp.dge
  ),
  file = "ov_GSE32062_dge_result.qs",
  nthreads = 4L
)

cli::cli_alert_success("All done!")

# PID='3258191'
data_path = "/home/data/sigbridger/GSEA/ov"
dge_res = qs::qread(
  file.path(data_path, "ov_GSE32062_dge_result.qs"),
  nthreads = 4L
)
dge_res$scab = scab.dge
qs::qsave(
  dge_res,
  file.path(data_path, "ov_GSE32062_dge_result.qs"),
  nthreads = 4L
)
