setwd(file.path(
  usethis::proj_path(),
  "paper_plot/Figure1A/GSEA"
))

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)

# TCGA-LUAD survival cohort: Positive cells versus non-Positive cells.
cohort <- "survival_lung_TCGA_LUAD_merged_seurat"
methods <- c("scissor", "scAB", "scPAS", "scPP", "DEGAS", "LP_SGL")
focal_pathways <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_G2M_CHECKPOINT"
)

fgsea_results <- qs::qread(
  "../../../5_GSEA/test_lung/gsea_res_lung.qs",
  nthreads = 4L
)
degs <- qs::qread(
  "../../../5_GSEA/test_lung/degs_lung.qs",
  nthreads = 4L
)
hallmarks <- data.table::fread("../../../5_GSEA/geneset_hallmark.csv") |>
  distinct(gs_name, gene_symbol)

output_dir <- "mini/luad_survival"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Positron does not support fgsea's default forked workers.
BiocParallel::register(BiocParallel::SerialParam())

run_hallmark_gsea <- function(gene_rank) {
  GSEA(
    geneList = sort(gene_rank, decreasing = TRUE),
    TERM2GENE = hallmarks |> select(gs_name, gene_symbol),
    minGSSize = 15,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 0,
    seed = FALSE,
    verbose = FALSE
  )
}

save_dotplot <- function(gsea_result, method) {
  result_table <- as.data.frame(gsea_result@result) |>
    filter(!is.na(p.adjust), p.adjust < 0.05) |>
    arrange(p.adjust)

  if (nrow(result_table) == 0L) {
    return(invisible(NULL))
  }

  significant_result <- gsea_result
  significant_result@result <- result_table

  p <- dotplot(significant_result, showCategory = 10, split = ".sign") +
    facet_grid(. ~ .sign) +
    scale_color_gradient(low = "#386c9b", high = "#a02020") +
    labs(
      title = method,
      x = "Normalized enrichment score",
      y = NULL,
      color = "FDR"
    ) +
    theme_classic(base_size = 7) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      axis.text.y = element_text(size = 5),
      legend.position = "bottom",
      strip.background = element_blank()
    )

  ggsave(
    file.path(output_dir, paste0("luad_", method, "_hallmark_dotplot.png")),
    p,
    width = 5.2,
    height = 2.7,
    dpi = 400
  )
}

save_running_score <- function(gsea_result, method, pathway) {
  if (!pathway %in% gsea_result@result$ID) {
    return(invisible(NULL))
  }

  pathway_label <- sub("HALLMARK_", "", pathway) |>
    gsub("_", " ", x = _) |>
    tools::toTitleCase()

  p <- gseaplot2(
    gsea_result,
    geneSetID = pathway,
    title = paste(method, "—", pathway_label),
    base_size = 8,
    subplots = 1:3,
    pvalue_table = TRUE
  )

  filename_pathway <- sub("HALLMARK_", "", pathway) |>
    tolower()
  ggsave(
    file.path(
      output_dir,
      paste0("luad_", method, "_", filename_pathway, "_running_score.png")
    ),
    p,
    width = 4.6,
    height = 3.2,
    dpi = 400
  )
}

set.seed(4729)
gsea_by_method <- lapply(methods, function(method) {
  gene_rank <- degs[[cohort]][[method]]
  if (is.null(gene_rank)) {
    return(NULL)
  }

  gsea_result <- run_hallmark_gsea(gene_rank)
  if (is.null(gsea_result) || nrow(gsea_result@result) == 0L) {
    return(NULL)
  }

  write.csv(
    as.data.frame(gsea_result@result) |> arrange(p.adjust),
    file.path(
      output_dir,
      paste0("luad_", method, "_hallmark_gsea_results.csv")
    ),
    row.names = FALSE
  )
  save_dotplot(gsea_result, method)
  lapply(focal_pathways, function(pathway) {
    save_running_score(gsea_result, method, pathway)
  })
  gsea_result
})
names(gsea_by_method) <- methods

gsea_by_method <- gsea_by_method[!vapply(gsea_by_method, is.null, logical(1))]

format_pathway <- function(pathway) {
  tools::toTitleCase(gsub("_", " ", sub("HALLMARK_", "", pathway)))
}

# Compact 2 × 3 panels: one running enrichment-score trace per method.
make_running_panel <- function(pathway) {
  plots <- lapply(names(gsea_by_method), function(method) {
    gseaplot2(
      gsea_by_method[[method]],
      geneSetID = pathway,
      title = method,
      base_size = 7,
      subplots = 1
    ) +
      theme(
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6)
      )
  })

  pathway_label <- format_pathway(pathway)

  panel <- patchwork::wrap_plots(plots, ncol = 3) +
    patchwork::plot_annotation(title = pathway_label) &
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

  ggsave(
    file.path(
      output_dir,
      paste0(
        "luad_all_methods_",
        tolower(sub("HALLMARK_", "", pathway)),
        "_panel.png"
      )
    ),
    panel,
    width = 9,
    height = 4.6,
    dpi = 400
  )
}

invisible(lapply(focal_pathways, make_running_panel))

# NES is encoded by tile color; each tile reports its FDR-adjusted p-value.
focal_summary <- do.call(
  rbind,
  lapply(names(gsea_by_method), function(method) {
    result <- as.data.frame(gsea_by_method[[method]]@result)
    result <- result[result$ID %in% focal_pathways, c("ID", "NES", "p.adjust")]
    result$method <- method
    result
  })
) |>
  rename(pathway = ID, fdr = p.adjust) |>
  mutate(
    method = factor(method, levels = methods),
    pathway = factor(pathway, levels = rev(focal_pathways)),
    fdr_label = case_when(
      is.na(fdr) ~ "NA",
      fdr < 0.001 ~ "<0.001",
      TRUE ~ formatC(fdr, format = "f", digits = 3)
    ),
    nes_label = sprintf("NES %.2f", NES)
  )

p_heatmap <- ggplot(focal_summary, aes(x = method, y = pathway, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = paste(nes_label, "FDR", fdr_label, sep = "\n")),
    size = 2.7
  ) +
  scale_fill_gradient2(
    low = "#386c9b",
    mid = "white",
    high = "#a02020",
    midpoint = 0,
    name = "NES"
  ) +
  scale_y_discrete(labels = format_pathway) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

ggsave(
  file.path(output_dir, "luad_all_methods_focal_pathway_nes_fdr_heatmap.png"),
  p_heatmap,
  width = 5.8,
  height = 2.8,
  dpi = 400
)

write.csv(
  focal_summary,
  file.path(output_dir, "luad_all_methods_focal_pathway_nes_fdr.csv"),
  row.names = FALSE
)
