library(data.table)

setwd(file.path(usethis::proj_path(), "5_GSEA/all_viz/combined_table"))

le <- fread(file.path(usethis::proj_path(), "5_GSEA/all_viz/leading_edge_stats.csv"))
gsea <- fread(file.path(usethis::proj_path(), "5_GSEA/all_viz/gsea_res_all_stats.csv"))

keys <- c("screen_method", "dataset", "tumor", "pathway")

res <- merge(
  le, gsea,
  by = keys, all = TRUE,
  suffixes = c(".le", ".gsea")
)

# leadingEdge 列可能来自任一表，合并为一列
res[, leadingEdge := fifelse(is.na(leadingEdge.le), leadingEdge.gsea, leadingEdge.le)]
res[, c("leadingEdge.le", "leadingEdge.gsea") := NULL]

setcolorder(res, c(keys, setdiff(names(res), keys)))

fwrite(res, "leading_edge_all_stats.csv")
