library(Seurat)

setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_binary/brca/TNBC")
data_path = "/home/data/sigbridger/benchmark_binary/brca/TNBC"


GSE162228_scissor_result = qs::qread(
  file.path(data_path, "binary_GSE162228_tnbc_scissor_result.qs"),
  nthreads = 8L
)


GSE162228_scissor_pos <- colnames(GSE162228_scissor_result$scRNA_data)[
  which(GSE162228_scissor_result$scRNA_data$scissor == "Positive")
]
# GSE162228_scissor_neg <- colnames(GSE162228_scissor_result$scRNA_data)[
#     which(GSE162228_scissor_result$scRNA_data$scissor == "Negative")
# ]

GSE162228_scab_result = qs::qread(
  file.path(data_path, "binary_GSE162228_tnbc_scab_result.qs"),
  nthreads = 8L
)

GSE162228_scab_pos <- colnames(GSE162228_scab_result$scRNA_data)[
  which(GSE162228_scab_result$scRNA_data$scAB == "Positive")
]
# GSE162228_scab_neg <- colnames(GSE162228_scab_result$scRNA_data)[
#     which(GSE162228_scab_result$scRNA_data$scAB == "Other")
# ]

GSE162228_scpas_result = qs::qread(
  file.path(data_path, "binary_GSE162228_tnbc_scpas_result.qs"),
  nthreads = 8L
)

GSE162228_scpas_pos <- colnames(GSE162228_scpas_result$scRNA_data)[
  which(GSE162228_scpas_result$scRNA_data$scPAS == "Positive")
]
# GSE162228_scpas_neg <- colnames(GSE162228_scpas_result$scRNA_data)[
#     which(GSE162228_scpas_result$scRNA_data$scPAS == "Negative")
# ]

# ! scPP is suitable for this dataset
GSE162228_scpp_result = qs::qread(
  file.path(data_path, "binary_GSE162228_tnbc_scpp_result.qs"),
  nthreads = 8L
)

# GSE162228_scpp_pos <- colnames(GSE162228_scpp_result$scRNA_data)[
#     which(GSE162228_scpp_result$scRNA_data$scPP == "Positive")
# ]
# GSE162228_scpp_neg <- colnames(GSE162228_scpp_result$scRNA_data)[
#     which(GSE162228_scpp_result$scRNA_data$scPP == "Negative")
# ]

GSE162228_degas_result = qs::qread(
  file.path(data_path, "binary_GSE162228_tnbc_degas_result.qs"),
  nthreads = 8L
)

GSE162228_lp_sgl_result = qs::qread(
  file.path(data_path, "binary_GSE162228_tnbc_lp_sgl_result.qs"),
  nthreads = 8L
)

# ? ---- Sample info ----

all_cells = colnames(GSE162228_scissor_result$scRNA_data)
# > all_cells |> length()
# [1] 31872
seurat_tnbc_tumor = readRDS(
  "/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_TNBCTum.rds"
)

tnbc_tumor_cells = rownames(seurat_tnbc_tumor@meta.data)

# > nrow(seurat_tnbc_tumor@meta.data)
# [1] 30353

# ? ----- HER2 tumor cell Venn diagram ----

pos_venn = list(
  scissor = GSE162228_scissor_pos,
  scab = GSE162228_scab_pos,
  scpas = GSE162228_scpas_pos,
  # scpp = GSE162228_scpp_pos,
  tnbc_tumor = tnbc_tumor_cells,
  tnbc_cells = all_cells
)
set.seed(123)

venn_plot = ggVennDiagram::ggVennDiagram(
  x = pos_venn,
  category.names = c(
    "Scissor",
    "scAB",
    "scPAS",
    # "scPP",
    "TNBC tumor cells",
    "TNBC cells"
  ),
  set_color = c(
    "red",
    "blue",
    "#37ae00ff",
    # "#9c8200ff",
    "purple",
    "cyan"
  )
) +
  ggplot2::scale_fill_gradient(low = "white", high = "#ffb6b6ff")

# ? ----

seurat_merged = SigBridgeR::MergeResult(
  GSE162228_scab_result,
  GSE162228_scpas_result,
  GSE162228_scpp_result,
  GSE162228_scissor_result,
  GSE162228_lp_sgl_result,
  GSE162228_degas_result
)

qs::qsave(
  seurat_merged,
  file.path(data_path, "GSE162228_tnbc_merged_seurat.qs"),
  nthreads = 8L
)

seurat_merged <- qs::qread(
  file.path(data_path, "GSE162228_tnbc_merged_seurat.qs"),
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
    "LP_SGL"
    # ,    "PIPET"
  ),
  show_plot = TRUE,
  n_intersections = 50
)

upset <- upset$plot + cowplot::theme_cowplot()

ggplot2::ggsave("GSE162228_tnbc_Upset.png", upset, width = 11, height = 11)
