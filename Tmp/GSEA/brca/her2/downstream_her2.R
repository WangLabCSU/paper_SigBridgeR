# setwd(.rs.api.getActiveDocumentContext()$path |> dirname())
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/GSEA/brca/her2")
save_path = "/home/data/sigbridger/GSEA/brca/her2"
irgsea_score = qs::qread(
  "/home/data/sigbridger/GSEA/brca/her2/her2_irGSEA_score.qs",
  nthreads = 8L
)

# ! BULK- GSE42568
# ! SC- GSE161529 - her2
# ! phenotype - survival

her2_GSE42568_merged <- qs::qread(
  "/home/data/sigbridger/benchmark_data/brca/HER2/GSE42568_her2_merged_seurat.qs"
)

labeled_irgsea = SigBridgeR::MergeResult(irgsea_score, her2_GSE42568_merged)

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
  method = c("AUCell", "UCell", "singscore", "ssgsea")
)

cli::cli_alert_success("Done!")

scpas.dge <- irGSEA::irGSEA.integrate(
  object = labeled_irgsea,
  group.by = "scPAS",
  metadata = NULL,
  col.name = NULL,
  method = c("AUCell", "UCell", "singscore", "ssgsea")
)

cli::cli_alert_success("Done!")

scab.dge <- irGSEA::irGSEA.integrate(
  object = labeled_irgsea,
  group.by = "scAB",
  metadata = NULL,
  col.name = NULL,
  method = c("AUCell", "UCell", "singscore", "ssgsea")
)

cli::cli_alert_success("Done!")

scpp.dge <- irGSEA::irGSEA.integrate(
  object = labeled_irgsea,
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
  file = file.path(save_path, "her2_GSE42568_dge_result.qs"),
  nthreads = 4L
)
# PID='1934790'

dge_res = qs::qread(file.path(save_path, "her2_GSE42568_dge_result.qs"))
dge_res$scab = scab.dge
qs::qsave(
  dge_res,
  file.path(save_path, "her2_GSE42568_dge_result.qs"),
  nthreads = 4L
)
