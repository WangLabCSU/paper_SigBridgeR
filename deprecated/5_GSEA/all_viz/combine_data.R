# ! run this after `mean_score.R`
# ! all tissue
library(dplyr)
library(data.table)

setwd(file.path(usethis::proj_path(), "5_GSEA/all_viz"))

data_dir <- "/home/data/sigbridger/GSEA/mean_score_gsea"
files <- list.files(
  data_dir,
  full.names = TRUE,
  include.dirs = FALSE
)
names(files) <- basename(files) %>%
  tools::file_path_sans_ext() %>%
  gsub("_mean.*", "", .)

gsea_score_loaded <- purrr::imap(
  files,
  ~ {
    cli::cli_alert_info("{.file {.y}} loaded")
    loaded <- data.table::fread(.x)
    loaded[, data_name := .y]
  }
)
gsea_score_combinedd <- dplyr::bind_rows(gsea_score_loaded)
data.table::set(
  gsea_score_combinedd,
  j = c("phenotype"),
  value = gsub("_.*", "", gsea_score_combinedd[, data_name])
)
data.table::set(
  gsea_score_combinedd,
  j = c("sc"),
  value = gsub(
    ".*(her2|tnbc|ov|lung).*",
    "\\1",
    gsea_score_combinedd[, data_name],
    ignore.case = TRUE
  )
)
data.table::set(
  gsea_score_combinedd,
  j = c("bulk"),
  value = gsub("^([^_]*_){2}", "", gsea_score_combinedd[, data_name])
)

data.table::fwrite(
  gsea_score_combinedd,
  file.path(data_dir, "gsea_score_combined.csv")
)
