# setwd(.rs.api.getActiveDocumentContext()$path |> dirname())
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/GSEA/binary_brca/binary_tum_tnbc"
)


# library(irGSEA)
source("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/GSEA/irGSEA_bubble.R")

irgsea_score = qs::qread(
  "/home/data/sigbridger/GSEA/brca/tnbc/tnbc_irGSEA_score.qs",
  nthreads = 8L
)

dge = qs::qread(
  "/home/data/sigbridger/GSEA/binary_brca/tnbc/binary_tnbc_GSE42568_dge.qs",
  nthreads = 8L
)

filtered_dge <- purrr::map(
  dge,
  ~ purrr::map(
    .x,
    ~ {
      if ("p_val" %in% names(.x)) {
        .x %>%
          dplyr::filter(p_val <= 0.05 & p_val_adj <= 0.05)
      } else {
        .x %>%
          dplyr::filter(pvalue <= 0.05)
      }
    }
  )
)

# ! 添加换行符，GO太长了
truncated_dge <- purrr::map(
  filtered_dge,
  ~ purrr::map(
    .x,
    ~ {
      .x = dplyr::rename(.x, "Full_name" = "Name")
      need2truncate <- .x$Full_name

      .x$Name <- purrr::map_chr(need2truncate, function(GO) {
        if (stringr::str_length(GO) > 30) {
          parts <- strsplit(GO, "-")[[1]]
          middle_index <- ceiling(length(parts) / 2)
          paste0(
            paste(parts[1:middle_index], collapse = "-"),
            "-\n",
            paste(
              parts[(middle_index + 1):length(parts)],
              collapse = "-"
            )
          )
        } else {
          GO
        }
      })
      .x
    }
  )
)

# ! 太丑了
# htmaps <- purrr::map(
#     filtered_dge,
#     ~ {
#         lapply(names(.x), function(method) {
#             irGSEA.heatmap(
#                 object = .x,
#                 method = method,
#                 top = 50,
#                 show.geneset = NULL,
#                 heatmap.width = 30,
#                 significance.color = c("#CECECE", "#ff857eff"),
#                 cluster.color = ggsci::pal_igv(),
#                 direction.color = c("#8abdffff", "#ff857eff")
#             )
#         })
#     }
# )

bubbles = purrr::map(
  truncated_dge,
  ~ {
    lapply(names(.x), function(method) {
      if (nrow(.x[[method]]) < 2) {
        return(NULL)
      }
      irGSEA.bubble(
        object = .x,
        method = method,
        significance.color = c("#CECECE", "#ff857eff"),
        cluster.color = setNames(
          c(
            "#ff3333",
            "#386c9b",
            "#CECECE",
            "#CECECE"
          ),
          c("Positive", "Negative", "Neutral", "Other")
        ),
        direction.color = c("#8abdffff", "#ff857eff"),
        cluster_rows = FALSE,
        top = 50
      )
    })
  }
)

plot_dir = "binary_plot_GSE42568"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}
purrr::imap(bubbles, function(dataset_list, dataset_name) {
  if (is.null(dataset_list) || length(dataset_list) == 0) {
    return(NULL)
  }
  method_name <- c("AUCell", "UCell", "singscore", "ssgsea", "RRA")

  purrr::imap(dataset_list, function(plot_obj, i) {
    if (is.null(plot_obj)) {
      return(NULL)
    }

    filename <- paste0(
      "tnbc_",
      dataset_name,
      "_",
      method_name[i],
      "_bubble.pdf"
    )
    filepath <- file.path(plot_dir, filename)

    ggplot2::ggsave(
      filename = filepath,
      plot = dataset_list[[i]] +
        ggplot2::theme(plot.margin = ggplot2::margin(l = 50)),
      height = 15,
      width = 10,
      limitsize = FALSE
    )

    message("已保存: ", filepath)
    return(filepath)
  })
})
