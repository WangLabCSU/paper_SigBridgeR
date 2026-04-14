setwd(file.path(usethis::proj_path(), "3_negative_ctrl/combined"))


test100_data <- list.files(
  path = "..",
  pattern = "rep100.*\\.csv",
  recursive = TRUE,
  full.names = TRUE
)
bulk <- dirname(test100_data) %>% basename()
sc_type <- dirname(dirname(test100_data)) %>% basename()
pheno_type <- gsub(".*(binary|survival).*", "\\1", dirname(test100_data))

names(test100_data) <- paste0(pheno_type, "_", sc_type, "_", bulk)

test100_data_loaded <- purrr::imap(
  .x = test100_data,
  .f = function(path, name) {
    dt <- data.table::fread(path)
    dt[, data_name := name] %>%
      tidyr::separate(
        col = "data_name",
        remove = FALSE,
        into = c("pheno_type", "sc_type", "bulk")
      )
  }
)

test100_data_combined <- dplyr::bind_rows(test100_data_loaded)
data.table::fwrite(x = test100_data_combined, file = "combined_100_test.csv")
