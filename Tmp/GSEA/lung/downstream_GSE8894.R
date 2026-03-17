# setwd(.rs.api.getActiveDocumentContext()$path |> dirname())
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/GSEA/lung")
data_path = "/home/data/sigbridger/GSEA/lung/"
irgsea_score = qs::qread(
  "/home/data/sigbridger/GSEA/lung/luad_irGSEA_score.qs",
  nthreads = 8L
)

luad_GSE8894 = qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/GSE8894/GSE8894_luad_merged_seurat.qs",
  nthreads = 8L
)

labeled_irgsea = SigBridgeR::MergeResult(irgsea_score, luad_GSE8894)


# Wlicox test is perform to all enrichment score matrixes and gene sets
# with adjusted p value &lt; 0.05 are used to integrated through RRA.
# Among them, Gene sets with p value &lt; 0.05 are statistically
# significant and common differential in all gene sets enrichment analysis
# methods. All results are saved in a list.

# scissor.dge <- irGSEA::irGSEA.integrate(
#     object = labeled_irgsea,
#     group.by = "scissor",
#     metadata = NULL,
#     col.name = NULL,
#     method = c("AUCell", "UCell", "singscore", "ssgsea")
# )

cli::cli_alert_success("Skip scissor due to no Positive cells available!")

# scpas.dge <- irGSEA::irGSEA.integrate(
#     object = labeled_irgsea,
#     group.by = "scPAS",
#     metadata = NULL,
#     col.name = NULL,
#     method = c("AUCell", "UCell", "singscore", "ssgsea")
# )

cli::cli_alert_success("Skip scPAS due to no Positive cells available!")

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
    # scissor = scissor.dge,
    scab = scab.dge,
    # scpas = scpas.dge,
    scpp = scpp.dge
  ),
  file = "luad_GSE8894_dge_result.qs",
  nthreads = 4L
)

cli::cli_h1("Done!")
# PID='3615925'
dge_res = qs::qread(
  file.path(data_path, "luad_GSE8894_dge_result.qs"),
  nthreads = 4L
)
dge_res$scab = scab.dge
qs::qsave(
  dge_res,
  file.path(data_path, "luad_GSE8894_dge_result.qs"),
  nthreads = 4L
)
