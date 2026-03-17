# setwd(.rs.api.getActiveDocumentContext()$path |> dirname())
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/GSEA/lung")

irgsea_score = qs::qread(
  "/home/data/sigbridger/GSEA/lung/luad_irGSEA_score.qs",
  nthreads = 8L
)

luad_tcga = qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/TCGA-LUAD/tcga_luad_merged_seurat.qs",
  nthreads = 8L
)

labeled_irgsea = SigBridgeR::MergeResult(irgsea_score, luad_tcga)


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

# scpp.dge <- irGSEA::irGSEA.integrate(
#     object = labeled_irgsea,
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
  file = "luad_tcga_dge_result.qs",
  nthreads = 4L
)

cli::cli_h1("Done!")
# PID='3738473'
data_path = "/home/data/sigbridger/GSEA/lung"
dge_res = qs::qread(
  file.path(data_path, "luad_tcga_dge_result.qs"),
  nthreads = 4L
)
dge_res$scab = scab.dge
qs::qsave(
  dge_res,
  file.path(data_path, "luad_tcga_dge_result.qs"),
  nthreads = 4L
)
