# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(file.path(usethis::proj_path(), "3_neg_ctrl_ucell/neg_ad_test"))

library(dplyr)
set.seed(123L)
mirai::daemons(4L)


neg_score <- qs::qread("../ucell_neg_score/neg_ucell_score.qs", nthreads = 8L)

neg_score_binded <- lapply(neg_score, function(x) {
  base::do.call(base::cbind, args = x)
})

seurats <- list.files(
  c(
    "/home/data/sigbridger/benchmark_data",
    "/home/data/sigbridger/benchmark_binary"
  ),
  pattern = "^survival.*\\.qs|^binary.*\\.qs",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE,
)
names(seurats) <- basename(seurats) %>%
  tools::file_path_sans_ext() %>%
  stringr::str_remove("_merged_seurat") %>%
  tolower()


method_labels <- if (!file.exists("method_labels.qs")) {
  mirai::mirai_map(seurats, function(file) {
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


neg_score_nested <- purrr::imap(
  neg_score_binded,
  function(
    rep100_score_mat,
    name # lowercase
  ) {
    l <- mirai::mirai_map(
      seq_along(method_labels),
      function(i) {
        data_name <- names(method_labels)[i]
        if (expected_tissue == "luad") {
          expected_tissue <- "lung"
        }

        if (!grepl(expected_tissue, data_name)) {
          return(NULL)
        }

        data <- cbind(rep100_score_mat, method_labels[[i]]) |>
          tidyr::pivot_longer(
            cols = names(method_labels[[i]]),
            names_to = "screen_method",
            values_to = "neg_ucell_score"
          ) |>
          tidyr::unite(
            col = "cluster",
            "screen_method",
            "neg_ucell_score"
          ) |>
          dplyr::group_by(cluster) |>
          dplyr::summarise(dplyr::across(
            dplyr::everything(),
            \(x) mean(x, na.rm = TRUE)
          ))
        gc(verbose = FALSE)

        data
      },
      method_labels = method_labels,
      expected_tissue = name,
      rep100_score_mat = rep100_score_mat
    )[
      mirai::.progress
    ]

    names(l) <- names(method_labels)

    l
  },
  .progress = "Nesting"
)

source("merge_nested_list.R")
neg_score_nested_combined <- merge_nested_lists(neg_score_nested)

qs::qsave(
  neg_score_nested_combined,
  file = "neg_score_nested_combined.qs",
  nthreads = 8L
)

ad_test_ucell_score <- purrr::imap(
  neg_score_nested_combined,
  function(tbl, name) {
    cli::cli_alert_info("Processing {.val {name}}")

    tbl %>%
      tidyr::pivot_longer(
        cols = colnames(.)[-1],
        names_to = "neg_sample",
        values_to = "ucell_mean_score"
      ) %>%
      dplyr::mutate(screen_group = gsub("_[^_]*$", "", cluster)) %>%
      dplyr::group_by(screen_group) %>%
      dplyr::group_modify(
        ~ {
          score_list <- split(.x$ucell_mean_score, .x$cluster)

          # 确保至少有 2 个分组
          if (length(score_list) < 2) {
            return(tibble::tibble(
              n_samples = NA_real_,
              n_groups = NA_real_,
              AS_stat = NA_real_,
              T_AV_stat = NA_real_,
              p_value = NA_real_,
              test_name = NA_character_,
              n_ties = NA_real_,
              sig = NA_real_,
              warning = NA_real_,
              null_dist1 = NA_real_,
              null_dist2 = NA_real_,
              method = NA_character_,
              n_sim = NA_real_,
              message = "分组数不足"
            ))
          }

          # 执行 k-sample Anderson-Darling 检验
          test_res <- tryCatch(
            {
              kSamples::ad.test(score_list)
            },
            error = function(e) {
              return(NULL)
            }
          )

          if (is.null(test_res)) {
            return(tibble::tibble(
              ad_statistic = NA,
              p_value = NA,
              n_groups = length(score_list),
              message = "检验失败"
            ))
          }

          tibble::tibble(
            n_samples = test_res$N,
            n_groups = test_res$k,
            AS_stat = test_res$ad[2, 1], # version 2
            T_AV_stat = test_res$ad[2, 2],
            p_value = test_res$ad[2, 3],
            test_name = test_res$test.name,
            n_ties = test_res$n.ties,
            sig = test_res$sig,
            warning = test_res$warning,
            null_dist1 = test_res$null.dist1,
            null_dist2 = test_res$null.dist2,
            method = test_res$method,
            n_sim = test_res$Nsim,
            message = "OK"
          )
        }
      ) %>%
      dplyr::ungroup()
  }
)

all_test_tbl <- purrr::imap_dfr(
  ad_test_ucell_score,
  ~ {
    dplyr::mutate(
      .x,
      data_name = .y,
      significant = p_value < 0.05,
      sig_label = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  }
)

data.table::fwrite(all_test_tbl, "ad_test_ucell_score.csv")

mirai::daemons(0L)

cli::cli_h1("AD test done")
gc()
