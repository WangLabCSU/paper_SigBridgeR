data_dir <- "/home/data/sigbridger/benchmark_data/ad"

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

# * GSE28146
# label: geo_accession
# value: title

# * GSE39420
# geo_accession
# title

# * GSE109887
# label: geo_accession
# value: source_name_ch1/title

pheno_files_converted <- purrr::imap(pheno_files_loaded, function(d, name) {
  d$screen_label <- ifelse(grepl("Control", d$title, ignore.case = TRUE), 0, 1)
  d
})

purrr::walk2(pheno_files_converted, pheno_files, function(d, file) {
  qs::qsave(d, file, nthreads = 4L)
})
