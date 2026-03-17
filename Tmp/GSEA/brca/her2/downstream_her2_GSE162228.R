# setwd(.rs.api.getActiveDocumentContext()$path |> dirname())
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/GSEA/brca/her2")

options(future.globals.maxSize = 100000 * 1024^5)

save_path = "/home/data/sigbridger/GSEA/brca/her2"

irgsea_score = qs::qread(
  "/home/data/sigbridger/GSEA/brca/her2/her2_irGSEA_score.qs",
  nthreads = 8L
)

# ! BULK- GSE162228
# ! SC- GSE161529 - her2
# ! phenotype - survival

her2_GSE162228_merged <- qs::qread(
  "/home/data/sigbridger/benchmark_data/brca/HER2/GSE162228_her2_merged_seurat.qs",
  nthreads = 4L
)

labeled_irgsea = SigBridgeR::MergeResult(irgsea_score, her2_GSE162228_merged)

# Wlicox test is perform to all enrichment score matrixes and gene sets
# with adjusted p value &lt; 0.05 are used to integrated through RRA.
# Among them, Gene sets with p value &lt; 0.05 are statistically
# significant and common differential in all gene sets enrichment analysis
# methods. All results are saved in a list.
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

# scissor.dge <- irGSEA::irGSEA.integrate(
#     object = labeled_irgsea,
#     group.by = "scissor",
#     metadata = NULL,
#     col.name = NULL,
#     method = c("AUCell", "UCell", "singscore", "ssgsea")
# )

cli::cli_alert_success(
  "scissor has been skipped due to no positive cells available!"
)

# # ! scPP is not available
# scpp.dge <- irGSEA::irGSEA.integrate(
#     object = labeled_irgsea,
#     group.by = "scPP",
#     metadata = NULL,
#     col.name = NULL,
#     method = c("AUCell", "UCell", "singscore", "ssgsea")
# )
# cli::cli_alert_success("Done!")

qs::qsave(
  list(
    # scissor = scissor.dge,
    scab = scab.dge,
    scpas = scpas.dge
    # ,        scpp = scpp.dge
  ),
  file = file.path(save_path, "her2_GSE162228_dge_result.qs"),
  nthreads = 4L
)
# PID='3084235'

dge_res = qs::qread(
  "/home/data/sigbridger/GSEA/brca/her2/her2_GSE162228_dge_result.qs",
  nthreads = 4L
)
dge_res$scab = scab.dge
qs::qsave(
  dge_res,
  file = file.path(save_path, "her2_GSE162228_dge_result.qs"),
  nthreads = 4L
)
