# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  file.path(usethis::proj_path(), "4_pos_ctrl_ucell/esmat")
)

library(dplyr)
library(data.table)

# ? Load marker files
markers_dir <- paste0(
  "../../4_positive_ctrl/",
  c("brca", "luad", "ov", "ad", "fshd")
)
marker_files <- list.files(
  markers_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)
names(marker_files) <- basename(marker_files) %>% tools::file_path_sans_ext()

loaded_marker_files <- purrr::map(marker_files, data.table::fread)

loaded_marker_splited <- purrr::imap(
  loaded_marker_files,
  function(raw_markers, name) {
    if (grepl("survival", name)) {
      n_risk <- raw_markers[direction == "risk", .N]
      n_prot <- raw_markers[direction == "protective", .N]
    } else {
      # * binary
      n_risk <- raw_markers[logFC > 0, .N]
      n_prot <- raw_markers[logFC < 0, .N]
    }
    if (n_risk < 20 || n_prot < 20) {
      cli::cli_warn(
        "{.val {name}}: Risk markers:  {n_risk}, Protective markers: {n_prot}"
      )
      return(NULL)
    }

    # ? convert to list and use top 20 (P.val < 0.05 first, log FC second)
    # ? already sorted
    if (grepl("survival", name)) {
      top_risk <- raw_markers[direction == "risk"][order(-abs_logHR)][1:20] %>%
        dplyr::pull(gene)
      top_protective <- raw_markers[direction == "protective"][order(
        -abs_logHR
      )][
        1:20
      ] %>%
        dplyr::pull(gene)
    } else {
      top_risk <- raw_markers[logFC > 0][1:20]$gene
      top_protective <- raw_markers[logFC < 0][1:20]$gene
    }

    gene_list <- list(
      top_risk,
      top_protective
    )
    names(gene_list) <- paste0(name, c("_pos_ucell", "_neg_ucell"))
    gene_list
  },
  .progress = "splitting"
)

# ? Seurat object
data_dir <- "/home/data/sigbridger/benchmark_data"
seurat_path <- file.path(
  data_dir,
  c(
    "brca/seurat_tnbc.qs",
    "lung/luad_GSE123902_seurat.qs",
    "brca/seurat_her2.qs",
    "ov/hgsoc_GSE165897_seurat.qs",
    "ad/ad_GSE138852_seurat.qs",
    "fshd/fshd_GSE122873_seurat.qs"
  )
)
seurat <- purrr::map(
  seurat_path,
  qs::qread,
  nthreads = 8L,
  .progress = "Reading Seurat"
)
names(seurat) <- c(
  "tnbc",
  "lung",
  "her2",
  "ov",
  "ad",
  "fshd"
)

# ? To locate data
sc_bulk_map <- function(sc = character(), bulk = character(), seurat = list()) {
  if (sc == "her2") {
    cli::cli_alert_info("Seurat: {.field BRCA HER2}")
    return(seurat$her2)
  }
  if (sc == "tnbc") {
    cli::cli_alert_info("Seurat: {.field TNBC}")
    return(seurat$tnbc)
  }

  switch(
    bulk,
    "GSE3141" = ,
    "GSE8894" = ,
    "GSE31210" = ,
    "TCGA_LUAD" = {
      cli::cli_alert_info("Seurat: {.field LUAD}")
      seurat$lung
    },
    "GSE140082" = ,
    "GSE9891" = ,
    "GSE32062_GPL6480" = {
      cli::cli_alert_info("Seurat: {.field OV}")
      seurat$ov
    },
    "GSE109887" = ,
    "GSE28146" = ,
    "GSE39420" = {
      cli::cli_alert_info("Seurat: {.field AD}")
      seurat$ad
    },
    "GSE115650" = ,
    "GSE140261" = ,
    "GSE56787" = {
      cli::cli_alert_info("Seurat: {.field FSHD}")
      seurat$fshd
    },
    cli::cli_abort("No seurat found: bulk: {bulk}, sc: {sc}")
  )
}

# ? run UCell
BPPARAM <- BiocParallel::MulticoreParam(workers = 16L)

# * list of matrices, row:cell, col: index
ucell_res <- purrr::imap(
  loaded_marker_splited,
  function(gene_list, gene_list_name) {
    cli::cli_alert_info("Handling {.val {gene_list_name}}")

    sc <- gsub(
      ".*(her2|tnbc|ad|fshd|ov).*",
      "\\1",
      gene_list_name,
      ignore.case = TRUE
    )
    bulk <- gsub(".*(TCGA.*|GSE.*)", "\\1", gene_list_name, ignore.case = TRUE)

    seurat_i <- sc_bulk_map(sc = sc, bulk = bulk, seurat = seurat)

    gc(verbose = FALSE)

    UCell::ScoreSignatures_UCell(
      matrix = SeuratObject::LayerData(seurat_i, layer = "data"),
      features = gene_list,
      precalc.ranks = NULL,
      maxRank = 1500,
      w_neg = 1,
      name = "_UCell",
      chunk.size = 100,
      BPPARAM = BPPARAM,
      ties.method = "average",
      force.gc = FALSE
    )
  }
)

qs::qsave(ucell_res, file = "ucell_res.qs", nthreads = 16L)

# ucell_res <- qs::qread("ucell_res.qs")

cli::cli_h1("Done!")
gc()
