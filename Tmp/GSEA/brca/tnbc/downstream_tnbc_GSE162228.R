# setwd(.rs.api.getActiveDocumentContext()$path |> dirname())
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/GSEA/brca/tnbc")
save_path = "/home/data/sigbridger/GSEA/brca/tnbc"
options(future.globals.maxSize = 100000 * 1024^5)

irgsea_score = qs::qread(
  "/home/data/sigbridger/GSEA/brca/tnbc/tnbc_irGSEA_score.qs",
  nthreads = 8L
)

# ! BULK- GSE162228
# ! SC- GSE161529 - tnbc
# ! phenotype - survival

tnbc_GSE162228_merged <- qs::qread(
  "/home/data/sigbridger/benchmark_data/brca/TNBC/GSE162228_tnbc_merged_seurat.qs"
)

labeled_irgsea = SigBridgeR::MergeResult(irgsea_score, tnbc_GSE162228_merged)

# Wlicox test is perform to all enrichment score matrixes and gene sets
# with adjusted p value &lt; 0.05 are used to integrated through RRA.
# Among them, Gene sets with p value &lt; 0.05 are statistically
# significant and common differential in all gene sets enrichment analysis
# methods. All results are saved in a list.

# scpas.dge <- irGSEA::irGSEA.integrate(
#     object = labeled_irgsea,
#     group.by = "scPAS",
#     metadata = NULL,
#     col.name = NULL,
#     method = c("AUCell", "UCell", "singscore", "ssgsea")
# )

cli::cli_alert_success("skip scPAS due to all Neutral cells!")

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

cli::cli_alert_success("skip scissor due to all Neutral cells!")

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
    scab = scab.dge
    # ,        scpas = scpas.dge
    # ,scpp = scpp.dge
  ),
  file = "tnbc_GSE162228_dge_result.qs",
  nthreads = 4L
)
# PID='3169000'
dge_res = qs::qread(
  file.path(save_path, "tnbc_GSE162228_dge_result.qs"),
  nthreads = 4L
)
dge_res$scab = scab.dge
qs::qsave(
  dge_res,
  file = file.path(save_path, "tnbc_GSE162228_dge_result.qs"),
  nthreads = 4L
)
