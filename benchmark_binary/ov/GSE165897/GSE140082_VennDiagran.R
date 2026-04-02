# setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_data/ov/GSE165897")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(Seurat)

# !  ovarian cancer
data_path = '/home/data/sigbridger/benchmark_binary/ov/GSE165897'

GSE140082_scissor_result = qs::qread(
  file.path(data_path, "GSE165897_GSE140082_Scissor_result.qs"),
  nthreads = 4
)


GSE140082_scissor_pos <- colnames(GSE140082_scissor_result$scRNA_data)[
  which(GSE140082_scissor_result$scRNA_data$scissor == "Positive")
]
# GSE140082_scissor_neg <- colnames(GSE140082_scissor_result$scRNA_data)[
#     which(GSE140082_scissor_result$scRNA_data$scissor == "Negative")
# ]

GSE140082_scab_result = qs::qread(
  file.path(data_path, "GSE165897_GSE140082_scAB_result.qs"),
  nthreads = 4
)
GSE140082_scab_pos <- colnames(GSE140082_scab_result$scRNA_data)[
  which(GSE140082_scab_result$scRNA_data$scAB == "Positive")
]
# GSE140082_scab_neg <- colnames(GSE140082_scab_result$scRNA_data)[
#     which(GSE140082_scab_result$scRNA_data$scab == "Negative")
# ]

GSE140082_scpas_result = qs::qread(
  file.path(data_path, "GSE165897_GSE140082_scPAS_result.qs"),
  nthreads = 4
)

GSE140082_scpas_pos <- colnames(GSE140082_scpas_result$scRNA_data)[
  which(GSE140082_scpas_result$scRNA_data$scPAS == "Positive")
]
# GSE140082_scpas_neg <- colnames(GSE140082_scpas_result$scRNA_data)[
#     which(GSE140082_scpas_result$scRNA_data$scPAS == "Negative")
# ]

GSE140082_scpp_result = qs::qread(
  file.path(data_path, "GSE165897_GSE140082_scPP_result.qs"),
  nthreads = 4
)

# GSE140082_scpp_pos <- colnames(GSE140082_scpp_result$scRNA_data)[
#     which(GSE140082_scpp_result$scRNA_data$scPP == "Positive")
# ]
# GSE140082_scpp_neg <- colnames(GSE140082_scpp_result$scRNA_data)[
#     which(GSE140082_scpp_result$scRNA_data$scPP == "Negative")
# ]

GSE140082_pipet_result = qs::qread(
  file.path(data_path, "GSE165897_GSE140082_PIPET_result.qs"),
  nthreads = 4
)

GSE140082_degas_result = qs::qread(
  file.path(data_path, "GSE165897_GSE140082_DEGAS_result.qs"),
  nthreads = 4
)

GSE140082_lp_sgl_result = qs::qread(
  file.path(data_path, "GSE165897_GSE140082_LP_SGL_result.qs"),
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
  scissor = GSE140082_scissor_pos,
  scpas = GSE140082_scpas_pos,
  scab = GSE140082_scab_pos,
  # scpp = GSE140082_scpp_pos,
  eoc_cells = eoc_cells,
  all_cells = all_cells
)
set.seed(123)

venn_plot = ggVennDiagram::ggVennDiagram(
  x = pos_venn,
  category.names = c(
    "Scissor",
    "scPAS",
    "scAB",
    # "scPP",
    "Epithelial ovarian cancer cells",
    "All cells"
  ),
  set_color = c(
    "red",
    "#37ae00ff",
    "blue",
    #  "#9c8200ff",
    "purple",
    "cyan"
  )
) +
  ggplot2::scale_fill_gradient(low = "white", high = "#ffb6b6ff") +
  ggplot2::ggtitle("GSE140082  ovarian cancer cells Venn diagram")

ggplot2::ggsave("GSE140082_OV_Venn.png", venn_plot, width = 10, height = 9)

# ? ----

seurat_merged = SigBridgeR::MergeResult(
  GSE140082_scab_result,
  GSE140082_scpp_result,
  GSE140082_scpas_result,
  GSE140082_scissor_result,
  GSE140082_lp_sgl_result,
  GSE140082_degas_result,
  GSE140082_pipet_result
)

qs::qsave(
  seurat_merged,
  file.path(data_path, "GSE140082_ov_merged_seurat.qs"),
  nthreads = 4
)

seurat_merged <- qs::qread(
  file.path(data_path, "GSE140082_ov_merged_seurat.qs"),
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

ggplot2::ggsave("GSE140082_ov_Upset.png", upset, width = 11, height = 11)
