# ! Negative Comparison
# ! Survival
# ! LUAD
# ! GSE8894
# ! GSE123902

setwd(
  '/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_negative_compare/survival/luad/GSE8894'
)

score_path <- '/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_negative_compare/luad'
score_name <- paste0('luad_Sample_', 1:10, '_ssgsea_score')
score_files <- paste0(score_name, '.qs')

# ? read files
purrr::walk2(
  score_files,
  score_name,
  function(file, name) {
    assign(name, qs::qread(file.path(score_path, file)), envir = .GlobalEnv)
  },
  .progress = TRUE
)

# ? A function to calculate the mean score of each group
MeanScore = function(seurat_obj, es.mat, sample = "GSE", i = NULL) {
  library(tidyr)

  meta <- seurat_obj[[]] %>%
    tibble::rownames_to_column("cell")

  # 处理表达矩阵并合并元数据
  es.df <- Matrix::t(es.mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::left_join(meta, by = "cell")

  score_row_name = if (is.null(i)) {
    'Neg_ctrl_Sample'
  } else {
    paste0('Neg_ctrl_Sample_', i)
  }

  mean_scores_scissor <- if ("scissor" %in% colnames(es.df)) {
    es.df %>%
      dplyr::group_by(scissor) %>%
      dplyr::summarise(dplyr::across(
        score_row_name,
        mean,
        na.rm = TRUE
      )) %>%
      pivot_longer(
        cols = score_row_name,
        names_to = "cluster",
        values_to = "mean_score"
      ) %>%
      pivot_longer(
        cols = "scissor",
        names_to = "screen_method",
        values_to = "screen_ressult"
      )
  } else {
    NULL
  }
  mean_scores_scpas <- if ("scPAS" %in% colnames(es.df)) {
    es.df %>%
      dplyr::group_by(scPAS) %>%
      dplyr::summarise(dplyr::across(
        score_row_name,
        mean,
        na.rm = TRUE
      )) %>%
      pivot_longer(
        cols = score_row_name,
        names_to = "cluster",
        values_to = "mean_score"
      ) %>%
      pivot_longer(
        cols = "scPAS",
        names_to = "screen_method",
        values_to = "screen_ressult"
      )
  } else {
    NULL
  }
  mean_scores_scab <- if ("scAB" %in% colnames(es.df)) {
    es.df %>%
      dplyr::group_by(scAB) %>%
      dplyr::summarise(dplyr::across(
        score_row_name,
        mean,
        na.rm = TRUE
      )) %>%
      pivot_longer(
        cols = score_row_name,
        names_to = "cluster",
        values_to = "mean_score"
      ) %>%
      pivot_longer(
        cols = "scAB",
        names_to = "screen_method",
        values_to = "screen_ressult"
      )
  } else {
    NULL
  }
  mean_scores_scpp <- if ("scPP" %in% colnames(es.df)) {
    es.df %>%
      dplyr::group_by(scPP) %>%
      dplyr::summarise(dplyr::across(
        score_row_name,
        mean,
        na.rm = TRUE
      )) %>%
      pivot_longer(
        cols = score_row_name,
        names_to = "cluster",
        values_to = "mean_score"
      ) %>%
      pivot_longer(
        cols = "scPP",
        names_to = "screen_method",
        values_to = "screen_ressult"
      )
  } else {
    NULL
  }

  mean_scores <- rbind(
    mean_scores_scissor,
    mean_scores_scpas,
    mean_scores_scab,
    mean_scores_scpp
  ) %>%
    unite("screen", screen_method, screen_ressult, sep = "_") %>%
    dplyr::mutate("Sample" = sample)
}

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/GSE8894/GSE8894_luad_merged_seurat.qs",
  nthreads = 4L
)

# ? Calculate the mean score of each group
mean_scores <- purrr::imap(
  score_name,
  function(name, index) {
    MeanScore(
      seurat_obj = seurat,
      es.mat = get(name),
      sample = 'GSE8894',
      i = index
    )
  },
  .progress = TRUE
)

names(mean_scores) <- score_name
tbl_mean_scores <- dplyr::bind_rows(mean_scores)

# ? Convert the cluster name to factor to make the order of the cluster fixed
tbl_mean_scores$cluster <- factor(
  tbl_mean_scores$cluster,
  levels = unique(tbl_mean_scores$cluster)
)


bubble <- ggplot2::ggplot(
  tbl_mean_scores,
  ggplot2::aes(x = cluster, y = screen, fill = mean_score)
) +
  ggplot2::geom_point(size = 6, alpha = 0.9, shape = 21, color = "black") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    name = "Mean Score"
  ) +
  ggplot2::labs(
    title = "Mean ssGSEA Scores by Group",
    subtitle = "sc-GSE123902, bulk-GSE8894",
    y = "Mean ssGSEA Score",
    x = "Group"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

ggplot2::ggsave(
  filename = 'Neg_ctrl_luad_GSE8894_mean_score.png',
  bubble,
  dpi = 400,
  width = 6,
  height = 5
)
