convert_gene_symbol <- function(gene) {
  genes <- gene |>
    stringr::str_remove_all("\\..*$")
  options(
    IDConverter.datapath = system.file("extdata", package = "IDConverter")
  )
  gene_symbols <- IDConverter::convert_hm_genes(genes, type = "ensembl")

  duplicated_or_na <- is.na(gene_symbols) | duplicated(gene_symbols)

  result <- gene_symbols

  if (any(duplicated_or_na)) {
    cli::cli_alert_info(
      "Some genes have multiple gene symbols or no gene symbol: {.val {which(duplicated_or_na)}}"
    )
    cli::cli_alert_info(
      "They are kept as is (ENSEMBL IDs)"
    )
    result[duplicated_or_na] <- gene[duplicated_or_na]
  }

  result
}
