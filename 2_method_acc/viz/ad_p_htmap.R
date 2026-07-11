ad_p_htmap <- function(
  corr, # wide data: p-value matrix
  p_mat = matrix(),
  color = colorRampPalette(c("#DDDDDD", "#B2182B")),
  filename = "ad_p_htmap.png",
  width = 800,
  height = 800,
  res = 400,
  tolerance = 1e-300
) {
  # log10 spread + min-max normalize to 0-1 for balanced color differentiation
  corr_trans <- -log10(pmax(corr, tolerance))

  Cairo::CairoPNG(
    filename = filename,
    width = width,
    height = height,
    res = res
  )

  corrplot::corrplot(
    corr = corr_trans,
    method = "ellipse",
    type = "lower",
    col = color(100L),
    outline = "grey",
    diag = TRUE,
    tl.cex = 1.2,
    tl.col = "black",
    addCoef.col = "black",
    number.cex = 0.8,
    is.corr = FALSE
  )

  corrplot::corrplot(
    corr = corr_trans,
    method = "pie",
    type = "upper",
    col = color(100L),
    outline = "grey",
    diag = TRUE,
    tl.cex = 1,
    tl.col = "black",
    tl.pos = "d",
    p.mat = p_mat,
    sig.level = c(0.001, 0.01, 0.05),
    insig = "label_sig",
    pch.cex = 1.2,
    addCoef.col = "black",
    number.cex = 0.01,
    add = TRUE,
    is.corr = FALSE
  )

  corrplot::corrplot(
    corr = corr_trans,
    method = "number",
    type = "lower",
    col = color(100L),
    outline = "grey",
    diag = FALSE,
    tl.pos = "n",
    cl.pos = "n",
    addCoef.col = "black",
    number.cex = 0.9,
    add = TRUE,
    is.corr = FALSE
  )
  dev.off()
}
