setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(data.table)

esmat_root <- "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/normalized_esmat/survival/brca/her2"
esmat_files <- list.files(esmat_root, recursive = TRUE) %>%
  grep("ssGSEA_score_z\\.qs", ., value = TRUE)

bulks <- gsub(".*(GSE\\d+|TCGA_[A-Z]{4}).*", "\\1", esmat_files) %>% unique()

type_pheno <- "survival"


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

  if (inherits(data, "data.table")) {
    target_val <- data[cluster == target_group, .(z_ssgsea_score)][[1L]]

    rest_val <- data[cluster %chin% rest_groups, .(z_ssgsea_score)][[1L]]
  } else {
    target_val <- data$z_ssgsea_score[data$cluster == target_group]

    rest_val <- data$z_ssgsea_score[data$cluster %chin% rest_groups]
  }

  test <- stats::wilcox.test(
    target_val,
    rest_val,
    p.adjust.method = "BH",
    conf.int = TRUE,
    exact = FALSE
  )

  #   test <- pairwise.wilcox.test(
  #     x = data$z_ssgsea_score,
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
    # get all the methods when using this bulk
    screen_method_i <- gsub(
      ".*(scissor|scPAS|scAB|scPP).*",
      "\\1",
      files_of_bulk_i
    ) %>%
      unique()

    # conclude the result of each bulk
    combined <- data.frame()

    purrr::walk2(
      screen_method_i,
      files_of_bulk_i,
      function(method_i, file_i) {
        # long format esmat
        esmat_i <- qs::qread(file.path(esmat_root, file_i))

        group_type <- unique(esmat_i$cluster)

        group_n <- length(group_type)

        # wilcoxon test
        if (group_n == 2) {
          # scAB or only two groups in scissor, scPAS and scPP

          cli::cli_alert_info(
            "Handling {bulk_i} - {method_i} with {group_n} groups"
          )

          #  fall back if error occur
          #   fall_back <- data.frame(
          #     p.value = -1,
          #     parameter = -1,
          #     statistic = -1,
          #     null.value = -1,
          #     alternative = "Wrong",
          #     method = "Wrong",
          #     data.name = "Wrong",
          #     p.adjust.method = "Wrong",
          #     screen_method = method_i,
          #     bulk = bulk_i,
          #     type_pheno = type_pheno
          #   )

          results <- dplyr::bind_rows(
            one_vs_rest_wilcox(esmat_i, group_type[1]) %>%
              pw_to_long(
                method_i,
                bulk_i,
                type_pheno,
                group_type[1]
              ),
            one_vs_rest_wilcox(esmat_i, group_type[2]) %>%
              pw_to_long(
                method_i,
                bulk_i,
                type_pheno,
                group_type[2]
              )
          )
        } else if (group_n == 3) {
          cli::cli_alert_info(
            "Handling {bulk_i} - {method_i} with {group_n} groups"
          )

          # ? Three groups
          # scissor or scPAS or scPP in normal condition

          results <- dplyr::bind_rows(
            one_vs_rest_wilcox(esmat_i, paste0(method_i, "_Negative")) %>%
              pw_to_long(
                method_i = method_i,
                bulk_i = bulk_i,
                type_pheno = type_pheno,
                target_group = paste0(method_i, "_Negative")
              ),
            one_vs_rest_wilcox(esmat_i, paste0(method_i, "_Neutral")) %>%
              pw_to_long(
                method_i = method_i,
                bulk_i = bulk_i,
                type_pheno = type_pheno,
                target_group = paste0(method_i, "_Neutral")
              ),
            one_vs_rest_wilcox(esmat_i, paste0(method_i, "_Positive")) %>%
              pw_to_long(
                method_i = method_i,
                bulk_i = bulk_i,
                type_pheno = type_pheno,
                target_group = paste0(method_i, "_Positive")
              )
          )
        } else {
          cli::cli_warn(
            "{bulk_i} - {method_i} has {group_n} groups, which is not supported."
          )
          results <- data.frame()
        }

        combined <<- dplyr::bind_rows(combined, results)
      }
    )

    assign(paste0("diff_", bulk_i), combined, envir = .GlobalEnv)
  },
  .progress = TRUE
)

# ? save result
purrr::walk(
  bulks,
  function(bulk_i) {
    diff_i <- get(paste0("diff_", bulk_i))
    file_name <- paste0("z_diff_", bulk_i, "_survival.csv")
    fwrite(diff_i, file = file_name)
  }
)
