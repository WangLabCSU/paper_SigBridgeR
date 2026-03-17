# setwd(.rs.api.getActiveDocumentContext()$path |> dirname())
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/GSEA/binary_ov")

data_path = "/home/data/sigbridger/GSEA/ov"

irgsea_score = qs::qread(
  file.path(data_path, "ov_GSE140082_irGSEA_score.qs"),
  nthreads = 8L
)

# ! BULK- GSE140082
# ! SC- GSE165897 - ov
# ! phenotype - binary variable

# labeled_irgsea = SigBridgeR::MergeResult(irgsea_score, ov_GSE9891_merged)
labeled_irgsea = irgsea_score


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

output_dir = "/home/data/sigbridger/GSEA/binary_ov"

qs::qsave(
  list(
    scissor = scissor.dge,
    scab = scab.dge,
    scpas = scpas.dge,
    scpp = scpp.dge
  ),
  file = file.path(output_dir, "binary_ov_GSE140082_dge.qs"),
  nthreads = 4L
)
# PID='4171853'

cli::cli_h1("All Done!")

dge_res = qs::qread(file.path(output_dir, "binary_ov_GSE140082_dge.qs"))

dge_res$scab = scab.dge
qs::qsave(dge_res, file.path(output_dir, "binary_ov_GSE140082_dge.qs"))
