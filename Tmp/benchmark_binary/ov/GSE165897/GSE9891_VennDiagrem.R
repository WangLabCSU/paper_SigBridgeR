setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(Seurat)

# ! HGSOC
# ! bulk- GSE9891
# ! sc- GSE165897
# ! binary phenotype

data_path = "/home/data/sigbridger/benchmark_binary/ov/GSE165897"

GSE9891_scissor_result = qs::qread(
  file.path(data_path, "GSE165897_GSE9891_Scissor_result.qs"),
  nthreads = 4
)


GSE9891_scissor_pos <- colnames(GSE9891_scissor_result$scRNA_data)[
  which(GSE9891_scissor_result$scRNA_data$scissor == "Positive")
]
# GSE9891_scissor_neg <- colnames(GSE9891_scissor_result$scRNA_data)[
#     which(GSE9891_scissor_result$scRNA_data$scissor == "Negative")
# ]

GSE9891_scab_result = qs::qread(
  file.path(data_path, "GSE165897_GSE9891_scAB_result.qs"),
  nthreads = 4
)

GSE9891_scab_pos <- colnames(GSE9891_scab_result$scRNA_data)[
  which(GSE9891_scab_result$scRNA_data$scAB == "Positive")
]
# GSE9891_scab_neg <- colnames(GSE9891_scab_result$scRNA_data)[
#     which(GSE9891_scab_result$scRNA_data$scAB == "Negative")
# ]

GSE9891_scpas_result = qs::qread(
  file.path(data_path, "GSE165897_GSE9891_scPAS_result.qs"),
  nthreads = 4
)

GSE9891_scpas_pos <- colnames(GSE9891_scpas_result$scRNA_data)[
  which(GSE9891_scpas_result$scRNA_data$scPAS == "Positive")
]
# GSE9891_scpas_neg <- colnames(GSE9891_scpas_result$scRNA_data)[
#     which(GSE9891_scpas_result$scRNA_data$scPAS == "Negative")
# ]

GSE9891_scpp_result = qs::qread(
  file.path(data_path, "GSE165897_GSE9891_scPP_result.qs"),
  nthreads = 4
)

GSE9891_scpp_pos <- colnames(GSE9891_scpp_result$scRNA_data)[
  which(GSE9891_scpp_result$scRNA_data$scPP == "Positive")
]
# GSE9891_scpp_neg <- colnames(GSE9891_scpp_result$scRNA_data)[
#     which(GSE9891_scpp_result$scRNA_data$scPP == "Negative")
# ]

GSE9891_DEGAS_result = qs::qread(
  file.path(data_path, "GSE165897_GSE9891_DEGAS_result.qs"),
  nthreads = 4
)


GSE9891_LP_SGL_result = qs::qread(
  file.path(data_path, "GSE165897_GSE9891_LP_SGL_result.qs"),
  nthreads = 4
)

GSE9891_PIPET_result = qs::qread(
  file.path(data_path, "GSE165897_GSE9891_PIPET_result.qs"),
  nthreads = 4
)


# ? ---- Sample info ----

seurat_GSE165897 = qs::qread(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_data/ov/hgsoc_GSE165897_seurat.qs",
  nthreads = 4
)

all_cells = colnames(seurat_GSE165897)

# > all_cells |>  length()
# [1] 50865

eoc_cells = colnames(seurat_GSE165897)[which(
  seurat_GSE165897$cell_type == "EOC"
)]


# ? ----- TNBC tumor cell Venn diagram ----

pos_venn = list(
  scissor = GSE9891_scissor_pos,
  scab = GSE9891_scab_pos,
  scpas = GSE9891_scpas_pos,
  scpp = GSE9891_scpp_pos,
  eoc_cells = eoc_cells,
  all_cells = all_cells
)
set.seed(123)

venn_plot = ggVennDiagram::ggVennDiagram(
  x = pos_venn,
  category.names = c(
    "Scissor",
    "scAB",
    "scPAS",
    "scPP",
    "Epithelial ovarian cancer cells",
    "All cells"
  ),
  set_color = c("red", "blue", "#37ae00ff", "#9c8200ff", "purple", "cyan")
) +
  ggplot2::scale_fill_gradient(low = "white", high = "#ffb6b6ff") +
  ggplot2::ggtitle("GSE9891 HGSOC cells Venn diagram")

ggplot2::ggsave("GSE9891_HGSOC_Venn.png", venn_plot, width = 10, height = 9)

# ? ----

seurat_merged = SigBridgeR::MergeResult(
  GSE9891_scab_result,
  GSE9891_scpp_result,
  GSE9891_scpas_result,
  GSE9891_scissor_result,
  GSE9891_LP_SGL_result,
  GSE9891_PIPET_result,
  GSE9891_DEGAS_result
)

qs::qsave(
  seurat_merged,
  file.path(data_path, "GSE9891_hgsoc_merged_seurat.qs"),
  nthreads = 4
)

seurat_merged <- qs::qread(
  file.path(data_path, "GSE9891_hgsoc_merged_seurat.qs"),
  nthreads = 4L
)

upset <- SigBridgeR::ScreenUpset(
  seurat_merged,
  screen_type = c(
    "scissor",
    "scAB",
    "scPAS",
    "scPP",
    "DEGAS",
    "LP_SGL",
    "PIPET"
  ),
  show_plot = TRUE,
  n_intersections = 50
)

upset <- upset$plot + cowplot::theme_cowplot()

ggplot2::ggsave("GSE9891_ov_Upset.png", upset, width = 11, height = 11)
