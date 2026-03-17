# setwd(.rs.api.getActiveDocumentContext()$path |> dirname())
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/GSEA/binary_brca/binary_tum_tnbc"
)

data_path = "/home/data/sigbridger/GSEA/brca/tnbc"

irgsea_score = qs::qread(
  file.path(data_path, "tnbc_irGSEA_score.qs"),
  nthreads = 8L
)

# ! BULK- tcga
# ! SC- GSE161529 - tnbc
# ! phenotype - binary variable

tnbc_tcga_merged <- qs::qread(
  "/home/data/sigbridger/benchmark_binary/brca/TNBC/tcga_tnbc_merged_seurat.qs",
  nthreads = 8L
)


labeled_irgsea = SigBridgeR::MergeResult(irgsea_score, tnbc_tcga_merged)

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

cli::cli_alert_success("skip scissor due to no Positive cells!")

scpas.dge <- irGSEA::irGSEA.integrate(
  object = labeled_irgsea,
  group.by = "scPAS",
  metadata = NULL,
  col.name = NULL,
  method = c("AUCell", "UCell", "singscore", "ssgsea")
)

cli::cli_alert_success("skip scPAS due to no Positive cells!")

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
cli::cli_alert_success("Skip scPP due to no label provided!")

output_dir = "/home/data/sigbridger/GSEA/binary_brca/tnbc"

qs::qsave(
  list(
    scissor = scissor.dge,
    scab = scab.dge,
    scpas = scpas.dge,
    scpp = scpp.dge
  ),
  file = file.path(output_dir, "binary_tnbc_tcga_dge.qs"),
  nthreads = 4L
)
# PID='683167'

cli::cli_h1("All Done!")
