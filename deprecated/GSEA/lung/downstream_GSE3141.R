# setwd(.rs.api.getActiveDocumentContext()$path |> dirname())
setwd(file.path(usethis::proj_path(), "Tmp/GSEA/lung"))

data_path <- "/home/data/sigbridger/GSEA/lung"

irgsea_score <- qs::qread(
  "/home/data/sigbridger/GSEA/lung/luad_irGSEA_score.qs",
  nthreads = 8L
)

luad_GSE3141 <- qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/GSE3141/GSE3141_luad_merged_seurat.qs",
  nthreads = 8L
)

irgsea_score <- SigBridgeR::MergeResult(irgsea_score, luad_GSE3141)


screen_labels <- grepv(
  "^sc[a-zA-Z]+$|DEGAS$|LP_SGL$|PIPET$",
  colnames(irgsea_score[[]])
)

if (file.exists(file.path(data_path, "luad_GSE3141_dge_result.qs"))) {
  dge_res <- qs::qread(
    file.path(data_path, "luad_GSE3141_dge_result.qs"),
    nthreads = 4L
  )

  done_labels <- names(dge_res)
} else {
  dge_res <- list()

  done_labels <- character(0)
}

is_more_than_2_group <- function(col, seurat) {
  if (length(unique(seurat@meta.data[[col]])) > 2) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Wlicox test is perform to all enrichment score matrixes and gene sets
# with adjusted p value &lt; 0.05 are used to integrated through RRA.
# Among them, Gene sets with p value &lt; 0.05 are statistically
# significant and common differential in all gene sets enrichment analysis
# methods. All results are saved in a list.

if (
  !"scissor" %in% done_labels &&
    "scissor" %in% screen_labels &&
    is_more_than_2_group("scissor", irgsea_score)
) {
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

if (
  !"scpas" %in% done_labels &&
    "scPAS" %in% screen_labels &&
    is_more_than_2_group("scPAS", irgsea_score)
) {
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

if (
  !"scab" %in% done_labels &&
    "scAB" %in% screen_labels &&
    is_more_than_2_group("scAB", irgsea_score)
) {
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

if (
  !"scpp" %in% done_labels &&
    "scPP" %in% screen_labels &&
    is_more_than_2_group("scPP", irgsea_score)
) {
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

if (
  !"lp_sgl" %in% done_labels &&
    "LP_SGL" %in% screen_labels &&
    is_more_than_2_group("LP_SGL", irgsea_score)
) {
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

if (
  !"degas" %in% done_labels &&
    "DEGAS" %in% screen_labels &&
    is_more_than_2_group("DEGAS", irgsea_score)
) {
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

if (
  !"pipet" %in% done_labels &&
    "PIPET" %in% screen_labels &&
    is_more_than_2_group("PIPET", irgsea_score)
) {
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

qs::qsave(
  dge_res,
  file.path(data_path, "luad_GSE3141_dge_result.qs"),
  nthreads = 4L
)
