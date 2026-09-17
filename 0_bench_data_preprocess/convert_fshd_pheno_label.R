data_dir <- "/home/data/sigbridger/benchmark_data/fshd"

pheno_files <- list.files(data_dir, full.names = TRUE, pattern = "pheno")
names(pheno_files) <- basename(pheno_files) %>%
  tools::file_path_sans_ext() %>%
  stringr::str_remove_all(".*_")

pheno_files_loaded <- lapply(pheno_files, qs::qread)

bulk_files <- list.files(data_dir, full.names = TRUE, pattern = "bulkdata")
names(bulk_files) <- basename(bulk_files) %>%
  tools::file_path_sans_ext() %>%
  stringr::str_remove_all(".*_")

bulk_files_loaded <- lapply(bulk_files, qs::qread)
col_names_type <- lapply(bulk_files_loaded, function(x) colnames(x) %>% head())
str(col_names_type) # all GSM.*

# * GSE26852
# label: geo_accession
# value: title

# * GSE56787
# geo_accession
# value: characteristics_ch1_1

# * GSE115650
# label: geo_accession
# value: source_name_ch1

# * GSE140261
# label: geo_accession
# value: ch1_disease state/title

pheno_files_converted <- purrr::imap(pheno_files_loaded, function(d, name) {
  d$screen_label <- if (name == "GSE26852") {
    ifelse(grepl("normal", d$title, ignore.case = TRUE), 0, 1)
  } else if (name == "GSE56787") {
    ifelse(grepl("control", d$characteristics_ch1_1, ignore.case = TRUE), 0, 1)
  } else if (name == "GSE115650") {
    ifelse(grepl("FSHD", d$source_name_ch1, ignore.case = TRUE), 1, 0)
  } else if (name == "GSE140261") {
    ifelse(grepl("control", d$title, ignore.case = TRUE), 0, 1)
  } else {
    stop("Unknown dataset")
  }

  d
})

purrr::walk2(pheno_files_converted, pheno_files, function(d, file) {
  qs::qsave(d, file, nthreads = 4L)
  cli::cli_alert_success("Saved file: {.file {file}}")
})
