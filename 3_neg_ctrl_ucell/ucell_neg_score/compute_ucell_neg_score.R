# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(file.path(usethis::proj_path(), "3_neg_ctrl_ucell/ucell_neg_score"))

set.seed(123L)

random_markers <- list.files(
  path = "../..",
  pattern = "random20_markers_100rep\\.csv",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)
names(random_markers) <- basename(random_markers) %>%
  tools::file_path_sans_ext() %>%
  gsub("_.*", "", .)

loaded_random_markers <- lapply(random_markers, data.table::fread)
