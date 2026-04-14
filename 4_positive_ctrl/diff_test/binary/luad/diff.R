setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(data.table)

esmat_root <- "../../../esmat/binary/luad"
esmat_files <- list.files(esmat_root, recursive = TRUE) %>%
  grep("ssGSEA_score.*\\.qs", ., value = TRUE)

bulks <- gsub(".*(GSE\\d+|TCGA_[A-Z]{4}).*", "\\1", esmat_files) %>% unique()

type_pheno <- "binary"
ssgsea_types <- c("pos_ssGSEA", "neg_ssGSEA")

# ? function to convert pw test result to long format
pw_to_long <- function(pw, method_i, bulk_i, type_pheno, target_group) {
  broom::tidy(pw) %>%
    mutate(
      parameter = pw$parameter,
      p.adjust.method = pw$p.adjust.method,
      null.value = pw$null.value,
      screen_method = method_i,
      bulk = bulk_i,
      type_pheno = type_pheno,
      group = paste0(target_group, " vs Rest")
    )
}

# ? function for one vs rest wilcox test
one_vs_rest_wilcox <- function(
  data,
  target_group
) {
  rest_groups <- setdiff(unique(data$cluster), target_group)

  target_val <- data[cluster == target_group, .(ssgsea_score)][[1L]]

  rest_val <- data[cluster %chin% rest_groups, .(ssgsea_score)][[1L]]

  test <- stats::wilcox.test(
    target_val,
    rest_val,
    p.adjust.method = "BH",
    conf.int = TRUE,
    exact = FALSE
  )

  #   test <- pairwise.wilcox.test(
  #     x = data$ssgsea_score,
  #     g = ifelse(
  #       data$cluster == target_group,
  #       target_group,
  #       "Rest"
  #     ),
  #     p.adjust.method = "BH",
  #     exact = FALSE
  #   )

  test
}


# ? diff test
purrr::walk(
  bulks,
  function(bulk_i) {
    # get esmat files of bulk_i
    files_of_bulk_i <- grep(bulk_i, esmat_files, value = TRUE)

    es_df <- qs::qread(file.path(esmat_root, files_of_bulk_i))
    # get all the methods when using this bulk
    screen_method <- gsub(
      ".*(scissor|scPAS|scAB|scPP|SCIPAC|LP_SGL|DEGAS|PIPET).*",
      "\\1",
      grepv(
        "scissor|scPAS|scAB|scPP|SCIPAC|LP_SGL|DEGAS|PIPET",
        colnames(es_df)
      )
    ) %>%
      unique()

    # conclude the result of each bulk
    combined <- data.frame()

    ## 🔁 Outer loop over ssGSEA types (pos / neg)
    purrr::walk(
      ssgsea_types,
      function(ss_type) {
        purrr::walk(
          screen_method,
          function(method_i) {
            if (!any(grepl(ss_type, colnames(es_df)))) {
              cli::cli_warn("No {ss_type} for {method_i}, skip")
              return(NULL)
            }

            # Select columns for current ssGSEA type + method column
            es_mat <- es_df %>%
              dplyr::select(
                dplyr::contains(ss_type), # <-- now dynamic: "pos_ssGSEA" or "neg_ssGSEA"
                !!sym(method_i)
              ) %>%
              tidyr::pivot_longer(
                dplyr::contains(ss_type), # <-- same dynamic here
                names_to = "ssGSEA type",
                values_to = "ssgsea_score"
              ) %>%
              dplyr::mutate(
                cluster = !!sym(method_i)
              ) %>%
              data.table::as.data.table()

            group_n <- length(unique(es_mat$cluster))
            group_type <- unique(es_mat$cluster)

            # wilcoxon test
            if (group_n == 2) {
              cli::cli_alert_info(
                "Handling {bulk_i} - {method_i} ({ss_type}) with {group_n} groups"
              )

              results <- dplyr::bind_rows(
                one_vs_rest_wilcox(es_mat, group_type[1]) %>%
                  pw_to_long(
                    method_i,
                    bulk_i,
                    type_pheno,
                    group_type[1]
                  ),
                one_vs_rest_wilcox(es_mat, group_type[2]) %>%
                  pw_to_long(
                    method_i,
                    bulk_i,
                    type_pheno,
                    group_type[2]
                  )
              ) %>%
                dplyr::mutate(
                  `ssGSEA type` = ss_type # <-- use current ss_type, not hardcoded "pos_ssGSEA"
                )
            } else if (group_n == 3) {
              cli::cli_alert_info(
                "Handling {bulk_i} - {method_i} ({ss_type}) with {group_n} groups"
              )

              results <- dplyr::bind_rows(
                one_vs_rest_wilcox(es_mat, group_type[1]) %>%
                  pw_to_long(
                    method_i = method_i,
                    bulk_i = bulk_i,
                    type_pheno = type_pheno,
                    target_group = group_type[1]
                  ),
                one_vs_rest_wilcox(es_mat, group_type[2]) %>%
                  pw_to_long(
                    method_i = method_i,
                    bulk_i = bulk_i,
                    type_pheno = type_pheno,
                    target_group = group_type[2]
                  ),
                one_vs_rest_wilcox(es_mat, group_type[3]) %>%
                  pw_to_long(
                    method_i = method_i,
                    bulk_i = bulk_i,
                    type_pheno = type_pheno,
                    target_group = group_type[3]
                  )
              ) %>%
                dplyr::mutate(
                  `ssGSEA type` = ss_type # <-- use current ss_type
                )
            } else {
              cli::cli_warn(
                "{bulk_i} - {method_i} ({ss_type}) has {group_n} groups, which is not supported."
              )
              results <- data.frame()
            }

            combined <<- dplyr::bind_rows(combined, results)
          }
        )
      }
    ) # end walk over ssgsea_types

    assign(paste0("diff_", bulk_i), combined, envir = .GlobalEnv)
  },
  .progress = TRUE
)

# ? save result
purrr::walk(
  bulks,
  function(bulk_i) {
    diff_i <- get(paste0("diff_", bulk_i))
    file_name <- paste0("diff_", bulk_i, "_binary.csv")
    fwrite(diff_i, file = file_name)
    cli::cli_alert_success("{file_name} saved")
  }
)
