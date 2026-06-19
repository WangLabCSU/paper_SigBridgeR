plot_gsea_table <- function(
  hallmarks = list(),
  fgseaRes = data.table::data.table(),
  deg = vector()
) {
  # up hallmarks
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n = 10), pathway]
  # down hallmarks
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n = 10), pathway]
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
    axisLabelStyle = NULL,
    render = NULL
  )
}
