library(Seurat)

setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_binary/brca/HER2")
data_path = "/home/data/sigbridger/benchmark_binary/brca/HER2"

tcga_scissor_result = qs::qread(
  file.path(data_path, "binary_her2_tcga_scissor_result.qs"),
  nthreads = 4
)


tcga_scissor_pos <- colnames(tcga_scissor_result$scRNA_data)[
  which(tcga_scissor_result$scRNA_data$scissor == "Positive")
]
# tcga_scissor_neg <- colnames(tcga_scissor_result$scRNA_data)[
#     which(tcga_scissor_result$scRNA_data$scissor == "Negative")
# ]

tcga_scab_result = qs::qread(
  file.path(data_path, "binary_her2_tcga_scab_result.qs"),
  nthreads = 4
)

tcga_scab_pos <- colnames(tcga_scab_result$scRNA_data)[
  which(tcga_scab_result$scRNA_data$scAB == "Positive")
]
# tcga_scab_neg <- colnames(tcga_scab_result$scRNA_data)[
#     which(tcga_scab_result$scRNA_data$scAB == "Other")
# ]

tcga_scpas_result = qs::qread(
  file.path(data_path, "binary_her2_tcga_scpas_result.qs"),
  nthreads = 4
)

tcga_scpas_pos <- colnames(tcga_scpas_result$scRNA_data)[
  which(tcga_scpas_result$scRNA_data$scPAS == "Positive")
]
# tcga_scpas_neg <- colnames(tcga_scpas_result$scRNA_data)[
#     which(tcga_scpas_result$scRNA_data$scPAS == "Negative")
# ]

tcga_scpp_result = qs::qread(
  file.path(data_path, "binary_her2_tcga_scpp_result.qs"),
  nthreads = 4
)

tcga_scpp_pos <- colnames(tcga_scpp_result$scRNA_data)[
  which(tcga_scpp_result$scRNA_data$scPP == "Positive")
]
# tcga_scpp_neg <- colnames(tcga_scpp_result$scRNA_data)[
#     which(tcga_scpp_result$scRNA_data$scPP == "Negative")
# ]

tcga_degas_result = qs::qread(
  file.path(data_path, "binary_her2_tcga_degas_result.qs"),
  nthreads = 4
)
tcga_lp_sgl_result = qs::qread(
  file.path(data_path, "binary_her2_tcga_lp_sgl_result.qs"),
  nthreads = 4
)
tcga_pipet_result = qs::qread(
  file.path(data_path, "binary_her2_tcga_pipet_result.qs"),
  nthreads = 4
)


# ? ---- Sample info ----

all_cells = colnames(tcga_scissor_result$scRNA_data)
# > all_cells |> length()
# [1] 31872
seurat_her2_tumor = readRDS(
  "/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_HER2Tum.rds"
)

her2_tumor_cells = rownames(seurat_her2_tumor@meta.data)

# > nrow(seurat_her2_tumor@meta.data)
# [1] 17849

# ? ----- HER2 tumor cell Venn diagram ----

pos_venn = list(
  scissor = tcga_scissor_pos,
  scab = tcga_scab_pos,
  scpas = tcga_scpas_pos,
  scpp = tcga_scpp_pos,
  her2_tumor = her2_tumor_cells,
  her2_cells = all_cells
)
set.seed(123)

venn_plot = ggVennDiagram::ggVennDiagram(
  x = pos_venn,
  category.names = c(
    "Scissor",
    "scAB",
    "scPAS",
    "scPP",
    "HER2 tumor cells",
    "HER2 cells"
  ),
  set_color = c(
    "red",
    "blue",
    "#37ae00ff",
    "#9c8200ff",
    "purple",
    "cyan"
  )
) +
  ggplot2::scale_fill_gradient(low = "white", high = "#ffb6b6ff")

# ? ----

seurat_merged = SigBridgeR::MergeResult(
  tcga_scab_result,
  tcga_scpas_result,
  tcga_scpp_result,
  tcga_scissor_result,
  tcga_lp_sgl_result,
  tcga_degas_result,
  tcga_pipet_result
)

qs::qsave(
  seurat_merged,
  file.path(data_path, "tcga_brca_merged_seurat.qs"),
  nthreads = 4L
)


seurat_merged <- qs::qread(
  file.path(data_path, "tcga_brca_merged_seurat.qs"),
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

ggplot2::ggsave("tcga_her2_Upset.png", upset, width = 11, height = 11)
