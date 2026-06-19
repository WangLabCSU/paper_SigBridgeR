plot_gsea_table <- function(
  hallmarks = list(),
  fgseaRes = data.table::data.table(),
  deg = vector()
) {
  ES <- pval <- pathway <- NULL

  # up hallmarks
  topPathwaysUp <- dtplyr::lazy_dt(fgseaRes) |>
    dplyr::filter(ES > 0 & !is.na(pval)) |>
    dplyr::arrange(pval) |>
    dplyr::slice(1:10) |>
    dplyr::pull(pathway)
  # down hallmarks
  topPathwaysDown <- dtplyr::lazy_dt(fgseaRes) |>
    dplyr::filter(ES < 0 & !is.na(pval)) |>
    dplyr::arrange(pval) |>
    dplyr::slice(1:10) |>
    dplyr::pull(pathway)

  topPathways <- c(topPathwaysUp, rev(topPathwaysDown)) # vector
  fgsea::plotGseaTable(
    pathways = hallmarks[topPathways],
    stats = deg,
    fgseaRes = fgseaRes,
    gseaParam = 0.5,
    colwidths = c(5, 3, 0.8, 1.2, 1.2),
    pathwayLabelStyle = NULL,
    headerLabelStyle = NULL,
    valueStyle = NULL,
    axisLabelStyle = list(size = 6),
    render = NULL
  )
  # ggplot2
}
