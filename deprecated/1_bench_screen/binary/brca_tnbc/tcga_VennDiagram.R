library(Seurat)

setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_binary/brca/TNBC")
data_path = "/home/data/sigbridger/benchmark_binary/brca/TNBC"


tcga_scissor_result = qs::qread(
  file.path(data_path, "binary_tnbc_tcga_scissor_result.qs"),
  nthreads = 4
)


tcga_scissor_pos <- colnames(tcga_scissor_result$scRNA_data)[
  which(tcga_scissor_result$scRNA_data$scissor == "Positive")
]
# tcga_scissor_neg <- colnames(tcga_scissor_result$scRNA_data)[
#     which(tcga_scissor_result$scRNA_data$scissor == "Negative")
# ]

tcga_scab_result = qs::qread(
  file.path(data_path, "binary_tnbc_tcga_scab_result.qs"),
  nthreads = 4
)

tcga_scab_pos <- colnames(tcga_scab_result$scRNA_data)[
  which(tcga_scab_result$scRNA_data$scAB == "Positive")
]
# tcga_scab_neg <- colnames(tcga_scab_result$scRNA_data)[
#     which(tcga_scab_result$scRNA_data$scAB == "Other")
# ]

tcga_scpas_result = qs::qread(
  file.path(data_path, "binary_tnbc_tcga_scpas_result.qs"),
  nthreads = 4
)

tcga_scpas_pos <- colnames(tcga_scpas_result$scRNA_data)[
  which(tcga_scpas_result$scRNA_data$scPAS == "Positive")
]
# tcga_scpas_neg <- colnames(tcga_scpas_result$scRNA_data)[
#     which(tcga_scpas_result$scRNA_data$scPAS == "Negative")
# ]

tcga_scpp_result = qs::qread(
  file.path(data_path, "binary_tnbc_tcga_scpp_result.qs"),
  nthreads = 4
)

# tcga_scpp_pos <- colnames(tcga_scpp_result$scRNA_data)[
#     which(tcga_scpp_result$scRNA_data$scPP == "Positive")
# ]
# tcga_scpp_neg <- colnames(tcga_scpp_result$scRNA_data)[
#     which(tcga_scpp_result$scRNA_data$scPP == "Negative")
# ]

# ? ---- Sample info ----

all_cells = colnames(tcga_scissor_result$scRNA_data)
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
  scissor = tcga_scissor_pos,
  scab = tcga_scab_pos,
  scpas = tcga_scpas_pos,
  # scpp = tcga_scpp_pos,
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
  tcga_scab_result,
  tcga_scpas_result,
  #   tcga_scpp_result,
  tcga_scissor_result,
  tcga_degas_result,
  tcga_lp_sgl_result,
  tcga_pipet_result
)

qs::qsave(
  seurat_merged,
  file.path(data_path, "tcga_tnbc_merged_seurat.qs"),
  nthreads = 4L
)

seurat_merged <- qs::qread(
  file.path(data_path, "tcga_tnbc_merged_seurat.qs"),
  nthreads = 4L
)

upset <- SigBridgeR::ScreenUpset(
  seurat_merged,
  screen_type = c(
    "scissor",
    "scAB",
    "scPAS",
    # "scPP",
    "DEGAS",
    "LP_SGL",
    "PIPET"
  ),
  show_plot = TRUE,
  n_intersections = 50
)

upset <- upset$plot + cowplot::theme_cowplot()

ggplot2::ggsave("tcga_tnbc_Upset.png", upset, width = 11, height = 11)
