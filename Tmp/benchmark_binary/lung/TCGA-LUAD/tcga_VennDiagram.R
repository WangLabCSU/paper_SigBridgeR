setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data_path = "/home/data/sigbridger/benchmark_binary/lung/TCGA-LUAD"

scissor_result = qs::qread(
  file.path(data_path, "luad_TCGA_LUAD_scissor.qs"),
  nthread = 4
)

tcga_scissor_pos <- colnames(scissor_result$scRNA_data)[
  which(scissor_result$scRNA_data$scissor == "Positive")
]

scpas_result = qs::qread(
  file.path(data_path, "luad_TCGA_LUAD_scpas.qs"),
  nthreads = 4
)

tcga_scpas_pos <- colnames(scpas_result$scRNA_data)[
  which(scpas_result$scRNA_data$scPAS == "Positive")
]

scab_result = qs::qread(
  file.path(data_path, "luad_TCGA_LUAD_scab.qs"),
  nthreads = 4
)

tcga_scab_pos <- colnames(scab_result$scRNA_data)[
  which(scab_result$scRNA_data$scAB == "Positive")
]

# ! scpp is not suitable for this dataset
# scpp_result = qs::qread("luad_TCGA_LUAD_scpp.qs", nthreads = 4)

# tcga_scpp_pos <- colnames(scpp_result$scRNA_data)[
#     which(scpp_result$scRNA_data$scPP == "Positive")
# ]

degas_result = qs::qread(
  file.path(data_path, "luad_TCGA_LUAD_degas.qs"),
  nthreads = 4
)

lp_sgl_result = qs::qread(
  file.path(data_path, "luad_TCGA_LUAD_lp_sgl.qs"),
  nthreads = 4
)

pipet_result = qs::qread(
  file.path(data_path, "luad_TCGA_LUAD_pipet.qs"),
  nthreads = 4
)


# ? ---- Sample info ----

all_cells = colnames(scissor_result$scRNA_data)

# > all_cells |>  length()
# [1] 32178

cnv_tumor_cells = colnames(scissor_result$scRNA_data)[which(
  scissor_result$scRNA_data$cnv_status == "tumor"
)]


# ? ----- TNBC tumor cell Venn diagram ----

pos_venn = list(
  scissor = tcga_scissor_pos,
  scpas = tcga_scpas_pos,
  scab = tcga_scab_pos,
  # scpp = tcga_scpp_pos,
  cnv_tumor_cells = cnv_tumor_cells,
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
    "CNV tumor cells",
    "All cells"
  ),
  set_color = c(
    "red",
    "#37ae00ff",
    "blue",
    # "#9c8200ff",
    "purple",
    "cyan"
  )
) +
  ggplot2::scale_fill_gradient(low = "white", high = "#ffb6b6ff") +
  ggplot2::ggtitle("TCGA LUAD cells Venn diagram")

ggplot2::ggsave("TCGA_LUAD_Venn.png", venn_plot, width = 10, height = 9)

# ? ----

seurat_merged = SigBridgeR::MergeResult(
  scpas_result,
  scissor_result,
  scab_result,
  # scpp_result,
  degas_result,
  lp_sgl_result,
  pipet_result
)

qs::qsave(
  seurat_merged,
  file.path(data_path, "tcga_luad_merged_seurat.qs"),
  nthreads = 4
)
seurat_merged <- qs::qread(
  file.path(data_path, "tcga_luad_merged_seurat.qs"),
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

ggplot2::ggsave("tcga_luad_Upset.png", upset, width = 11, height = 11)
