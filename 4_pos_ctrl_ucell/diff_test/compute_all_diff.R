setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(data.table)

esmat_files <- qs::qread("../esmat/ucell_res.qs", nthreads = 4L)

seurat_path <- c(
  "/home/data/sigbridger/benchmark_data/",
  "/home/data/sigbridger/benchmark_binary/"
)
seurat_files <- list.files(
  seurat_path,
  pattern = "^survival.*\\.qs|^binary.*\\.qs",
  full.names = TRUE,
  recursive = TRUE
)
names(seurat_files) <- basename(seurat_files) %>%
  tools::file_path_sans_ext() %>%
  gsub("_merged_seurat", "", .) %>%
  tolower()

mirai::daemons(4L)

# ? Read labels
method_labels <- if (!file.exists("method_labels.qs")) {
  mirai::mirai_map(seurat_files, function(file) {
    seurat_i <- qs::qread(file, nthreads = 8L)
    expected_methods <- c(
      "scissor",
      "scab",
      "scpas",
      "scpp",
      "degas",
      "lp_sgl",
      "pipet",
      "scipac"
    )
    meta <- seurat_i[[]] # data.frame
    col_names <- tolower(colnames(meta))

    meta[col_names %in% expected_methods]
  })[mirai::.progress]
  qs::qsave(method_labels, "method_labels.qs", nthreads = 8L)
} else {
  cli::cli_alert_info("Found existing {.val method_labels}")
  qs::qread("method_labels.qs", nthreads = 8L)
}

# ? "Other"/"Neutral"/"Negative" -> "Positive"
binarized_method_labels <- mirai::mirai_map(
  method_labels,
  function(df) {
    df %>%
      dplyr::mutate(dplyr::across(
        dplyr::everything(),
        ~ dplyr::case_when(
          . == "Positive" ~ "Positive",
          TRUE ~ "non_Positive"
        )
      ))
  },
  `%>%` = dplyr::`%>%`
)[mirai::.progress]


# ? To locate data
ucell_mat_map <- function(
  sc = character(),
  bulk = character(),
  pheno = character(),
  esmat_files = list()
) {
  bulk <- tolower(gsub("(tcga).*", "\\1", bulk, ignore.case = TRUE))
  names(esmat_files) <- tolower(names(esmat_files))

  chosen <- names(esmat_files) %>%
    grepv(pheno, .) %>%
    grepv(bulk, .)
  chosen <- if (!sc %in% c("lung", "ov")) {
    grepv(sc, chosen) # her2|tnbc
  } else {
    grepv("her2|tnbc", chosen, invert = TRUE)
  }
  if (length(chosen) == 0L || length(chosen) > 1L) {
    cli::cli_abort(c(
      "x" = "No unique match found or multiple matches found:\
       sc:{.val {sc}}, bulk:{.val {bulk}}, pheno:{.val {pheno}}",
      ">" = "chosen:{.val {chosen}}"
    ))
  }
  cli::cli_alert_success(
    "Detected sc: {.val {sc}}, bulk: {.val {bulk}}, pheno: {.val {pheno}}"
  )
  esmat_files[[chosen]]
}

# ? Find matched datasets
find_name <- function(chr = character) {
  sc <- gsub(
    ".*(her2|tnbc|lung|ov).*",
    "\\1",
    chr,
    ignore.case = TRUE
  )
  bulk <- tolower(gsub(
    ".*(TCGA.*|GSE.*)",
    "\\1",
    chr,
    ignore.case = TRUE
  ))
  pheno <- gsub(
    "(survival|binary).*",
    "\\1",
    chr,
    ignore.case = TRUE
  )
  list(
    sc = sc,
    bulk = bulk,
    pheno = pheno
  )
}


# ? Make data order the same
matched_data <- if (!file.exists("matched_data.qs")) {
  purrr::imap(
    binarized_method_labels,
    function(labels_of_seurat, names_of_seurat) {
      names_of_seurat <- tolower(names_of_seurat)

      l <- find_name(chr = names_of_seurat)
      sc <- l$sc
      bulk <- l$bulk
      pheno <- l$pheno

      cli::cli_alert_info(
        "Processing dataset: sc: {.val {sc}}, bulk: {.val {bulk}}, pheno: {.val {pheno}}"
      )

      ucell_mat <- ucell_mat_map(
        sc = sc,
        bulk = bulk,
        pheno = pheno,
        esmat_files = esmat_files
      )
      if (!all(row.names(ucell_mat) == row.names(labels_of_seurat))) {
        cli::cli_abort(
          "x" = "Row names of ucell_mat and labels_of_seurat are not the same: \
         {.val {names_of_seurat}}"
        )
      }

      cbind(ucell_mat, labels_of_seurat)
    }
  )
  qs::qsave(matched_data, "matched_data.qs", nthreads = 8L)
} else {
  cli::cli_alert_info("Found existing {.val matched_data}")
  qs::qread("matched_data.qs", nthreads = 8L)
}


all_diff <- mirai::mirai_map(
  seq_along(matched_data),
  function(i) {
    df <- matched_data[[i]]
    name_of_one_data <- names(matched_data)[i]

    expected_methods <- c(
      "scissor",
      "scab",
      "scpas",
      "scpp",
      "degas",
      "lp_sgl",
      "pipet",
      "scipac"
    )
    # First two columns are the two scores (pos_ucell, neg_ucell)
    score_names <- colnames(df)[1:2]
    # Find which columns correspond to methods
    existing_method_cols <- colnames(df)[
      tolower(colnames(df)) %in% expected_methods
    ]

    l <- find_name(chr = name_of_one_data)
    sc <- l$sc
    bulk <- l$bulk
    pheno <- l$pheno

    # Accumulate wilcox results for BH adjustment
    result_list <- vector(
      mode = "list",
      length = length(expected_methods) * length(score_names)
    )

    for (j in seq_along(existing_method_cols)) {
      method <- existing_method_cols[j]
      positive_idx <- df[[method]] == "Positive"
      non_positive_idx <- !positive_idx

      for (k in seq_along(score_names)) {
        score <- score_names[k]
        target_val <- df[[score]][positive_idx]
        rest_val <- df[[score]][non_positive_idx]

        wt <- tryCatch(
          stats::wilcox.test(
            target_val,
            rest_val,
            conf.int = TRUE,
            p.adjust.method = "BH",
            exact = FALSE
          ),
          error = function(e) NULL # only one group
        )

        mean_pos <- mean(target_val, na.rm = TRUE)
        mean_non <- mean(rest_val, na.rm = TRUE)

        res <- if (!is.null(wt)) {
          list(
            sc = sc,
            bulk = bulk,
            pheno = pheno,
            screen_method = method,
            score = score, # score name
            parameter = wt$parameter,
            null_value = wt$null.value,
            alternative = wt$alternative,
            method = wt$method,
            statistic = wt$statistic,
            p_value = wt$p.value,
            conf_int_low = wt$conf.int[1],
            conf_int_high = wt$conf.int[2],
            estimate = wt$estimate,
            mean_positive = mean_pos,
            mean_non_positive = mean_non,
            diff = mean_pos / mean_non
          )
        } else {
          list(
            sc = sc,
            bulk = bulk,
            pheno = pheno,
            screen_method = method,
            score = score, # score name
            parameter = NA_character_,
            null_value = NA_real_,
            alternative = NA_character_,
            method = NA_character_,
            statistic = NA_real_,
            p_value = NA_real_,
            conf_int_low = NA_real_,
            conf_int_high = NA_real_,
            estimate = NA_real_,
            mean_positive = mean_pos,
            mean_non_positive = mean_non,
            diff = mean_pos / mean_non
          )
        }

        result_list[[j * length(score_names) + k]] <- res
      }
    }

    dplyr::bind_rows(result_list)
  },
  matched_data = matched_data,
  find_name = find_name
)[mirai::.progress]

# Combine all datasets into one data.frame
all_diff_df <- dplyr::bind_rows(all_diff)

# Save result
data.table::fwrite(all_diff_df, "all_diff_df.csv")

cli::cli_h1("all_diff_df computed and saved")

mirai::daemons(0L)
gc()
