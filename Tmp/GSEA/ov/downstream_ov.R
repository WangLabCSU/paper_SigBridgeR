# setwd(.rs.api.getActiveDocumentContext()$path |> dirname())
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/GSEA/ov")

irgsea_score = qs::qread(
  "/home/data/sigbridger/GSEA/ov/ov_GSE140082_irGSEA_score.qs",
  nthreads = 8L
)

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

scpp.dge <- irGSEA::irGSEA.integrate(
  object = irgsea_score,
  group.by = "scPP",
  metadata = NULL,
  col.name = NULL,
  method = c("AUCell", "UCell", "singscore", "ssgsea")
)
cli::cli_alert_success("Done!")

qs::qsave(
  list(
    scissor = scissor.dge,
    scab = scab.dge,
    scpas = scpas.dge,
    scpp = scpp.dge
  ),
  file = "ov_GSE140082_dge_result.qs",
  nthreads = 4L
)
# PID='1750365'
data_path = "/home/data/sigbridger/GSEA/ov"
dge_res = qs::qread(
  file.path(data_path, "ov_GSE140082_dge_result.qs"),
  nthreads = 4L
)
dge_res$scab = scab.dge
qs::qsave(
  dge_res,
  file.path(data_path, "ov_GSE140082_dge_result.qs"),
  nthreads = 4L
)
