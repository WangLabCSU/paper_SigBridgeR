# setwd(.rs.api.getActiveDocumentContext()$path |> dirname())
setwd(file.path(usethis::proj_path(), "Tmp/GSEA/ov"))

# ? we just need score
irgsea_score <- qs::qread(
  "/home/data/sigbridger/GSEA/ov/ov_GSE140082_irGSEA_score.qs",
  nthreads = 8L
)

# # label
# my <- qs::qread(
#   "/home/data/sigbridger/benchmark_data/ov/GSE165897/GSE140082_ov_merged_seurat.qs",
#   nthreads = 4L
# )

ir_meta <- irgsea_score[[]]

# my_meta <- my[[]]

screen_labels <- grepv(
  "^sc[a-zA-Z]+$|DEGAS$|LP_SGL$|PIPET$",
  colnames(ir_meta)
)

# for (col in screen_labels) {
#   ir_meta[[col]] <- my_meta[col]
# }

# cli::cli_alert_info("{.field IrGSEA} Screen labels: {.val {screen_labels}}")

# irgsea_score@meta.data <- ir_meta

# qs::qsave(
#   irgsea_score,
#   "/home/data/sigbridger/GSEA/ov/ov_GSE140082_irGSEA_score.qs",
#   nthreads = 16L
# )

# Wlicox test is perform to all enrichment score matrixes and gene sets
# with adjusted p value &lt; 0.05 are used to integrated through RRA.
# Among them, Gene sets with p value &lt; 0.05 are statistically
# significant and common differential in all gene sets enrichment analysis
# methods. All results are saved in a list.
data_path <- "/home/data/sigbridger/GSEA/ov"


if (file.exists(file.path(data_path, "ov_GSE140082_dge_result.qs"))) {
  dge_res <- qs::qread(
    file.path(data_path, "ov_GSE140082_dge_result.qs"),
    nthreads = 4L
  )

  done_labels <- names(dge_res)
} else {
  dge_res <- list()

  done_labels <- character(0)
}

if (!"scissor" %in% done_labels && "scissor" %in% screen_labels) {
  scissor.dge <- irGSEA::irGSEA.integrate(
    object = irgsea_score,
    group.by = "scissor",
    metadata = NULL,
    col.name = NULL,
    method = c("AUCell", "UCell", "singscore", "ssgsea")
  )
  dge_res$scissor <- scissor.dge
  cli::cli_h2("Scissor Done!")
} else {
  cli::cli_alert_info("skip scissor")
}

if (!"scpas" %in% done_labels && "scPAS" %in% screen_labels) {
  scpas.dge <- irGSEA::irGSEA.integrate(
    object = irgsea_score,
    group.by = "scPAS",
    metadata = NULL,
    col.name = NULL,
    method = c("AUCell", "UCell", "singscore", "ssgsea")
  )
  dge_res$scpas <- scpas.dge
  cli::cli_h2("scPAS Done!")
} else {
  cli::cli_alert_info("skip scpas")
}

if (!"scab" %in% done_labels && "scAB" %in% screen_labels) {
  scab.dge <- irGSEA::irGSEA.integrate(
    object = irgsea_score,
    group.by = "scAB",
    metadata = NULL,
    col.name = NULL,
    method = c("AUCell", "UCell", "singscore", "ssgsea")
  )
  dge_res$scscabAB <- scab.dge
  cli::cli_h2("scAB Done!")
} else {
  cli::cli_alert_info("skip scab")
}

if (!"scpp" %in% done_labels && "scPP" %in% screen_labels) {
  scpp.dge <- irGSEA::irGSEA.integrate(
    object = irgsea_score,
    group.by = "scPP",
    metadata = NULL,
    col.name = NULL,
    method = c("AUCell", "UCell", "singscore", "ssgsea")
  )
  dge_res$scpp <- scpp.dge
  cli::cli_h2("scPP Done!")
} else {
  cli::cli_alert_info("skip scpp")
}

if (!"lp_sgl" %in% done_labels && "LP_SGL" %in% screen_labels) {
  lp_sgl.dge <- irGSEA::irGSEA.integrate(
    object = irgsea_score,
    group.by = "LP_SGL",
    metadata = NULL,
    col.name = NULL,
    method = c("AUCell", "UCell", "singscore", "ssgsea")
  )
  dge_res$lp_sgl <- lp_sgl.dge
  cli::cli_h2("LP_SGL Done!")
} else {
  cli::cli_alert_info("skip lp_sgl")
}

if (!"degas" %in% done_labels && "DEGAS" %in% screen_labels) {
  degas.dge <- irGSEA::irGSEA.integrate(
    object = irgsea_score,
    group.by = "DEGAS",
    metadata = NULL,
    col.name = NULL,
    method = c("AUCell", "UCell", "singscore", "ssgsea")
  )
  dge_res$degas <- degas.dge
  cli::cli_h2("DEGAS Done!")
} else {
  cli::cli_alert_info("skip degas")
}

if (!"pipet" %in% done_labels && "PIPET" %in% screen_labels) {
  pipet.dge <- irGSEA::irGSEA.integrate(
    object = irgsea_score,
    group.by = "PIPET",
    metadata = NULL,
    col.name = NULL,
    method = c("AUCell", "UCell", "singscore", "ssgsea")
  )
  dge_res$pipet <- pipet.dge
  cli::cli_h2("PIPET Done!")
} else {
  cli::cli_alert_info("skip pipet")
}

# ! store the survival label
qs::qsave(
  dge_res,
  file.path(data_path, "ov_GSE140082_dge_result.qs"),
  nthreads = 4L
)

res <- qs::qread(
  file.path(data_path, "ov_GSE140082_dge_result.qs"),
  nthreads = 4L
)
