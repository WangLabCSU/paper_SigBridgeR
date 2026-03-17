# setwd(.rs.api.getActiveDocumentContext()$path |> dirname())
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/GSEA/brca/tnbc")

options(future.globals.maxSize = 100000 * 1024^5)

save_path = "/home/data/sigbridger/GSEA/brca/tnbc"

irgsea_score = qs::qread(
  "/home/data/sigbridger/GSEA/brca/tnbc/tnbc_irGSEA_score.qs",
  nthreads = 8L
)

# ! BULK- tcga
# ! SC- GSE161529 - tnbc
# ! phenotype - survival

tnbc_tcga_merged <- qs::qread(
  "/home/data/sigbridger/benchmark_data/brca/TNBC/tcga_tnbc_merged_seurat.qs",
  nthreads = 4L
)

labeled_irgsea = SigBridgeR::MergeResult(irgsea_score, tnbc_tcga_merged)

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

cli::cli_alert_success("Done scPAS!")

scab.dge <- irGSEA::irGSEA.integrate(
  object = labeled_irgsea,
  group.by = "scAB",
  metadata = NULL,
  col.name = NULL,
  method = c("AUCell", "UCell", "singscore", "ssgsea")
)

cli::cli_alert_success("Done scAB!")

scissor.dge <- irGSEA::irGSEA.integrate(
  object = labeled_irgsea,
  group.by = "scissor",
  metadata = NULL,
  col.name = NULL,
  method = c("AUCell", "UCell", "singscore", "ssgsea")
)

cli::cli_alert_success(
  "Done scissor !"
)

scpp.dge <- irGSEA::irGSEA.integrate(
  object = labeled_irgsea,
  group.by = "scPP",
  metadata = NULL,
  col.name = NULL,
  method = c("AUCell", "UCell", "singscore", "ssgsea")
)
cli::cli_alert_success("Done scPP!")

qs::qsave(
  list(
    scissor = scissor.dge,
    scab = scab.dge,
    scpas = scpas.dge,
    scpp = scpp.dge
  ),
  file = file.path(save_path, "tnbc_tcga_dge_result.qs"),
  nthreads = 4L
)
# PID='575011'
