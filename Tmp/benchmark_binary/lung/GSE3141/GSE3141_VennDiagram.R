setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_data/lung/GSE3141")

library(Seurat)

# ! LUAD
# ! sc-GSE123902
# ! bulk-GSE3141

GSE3141_scissor_result = qs::qread(
  "luad_GSE3141_scissor.qs",
  nthreads = 4
)


GSE3141_scissor_pos <- colnames(GSE3141_scissor_result$scRNA_data)[
  which(GSE3141_scissor_result$scRNA_data$scissor == "Positive")
]
# GSE3141_scissor_neg <- colnames(GSE3141_scissor_result$scRNA_data)[
#     which(GSE3141_scissor_result$scRNA_data$scissor == "Negative")
# ]

GSE3141_scpas_result = qs::qread(
  "luad_GSE3141_scpas.qs",
  nthreads = 4
)

GSE3141_scpas_pos <- colnames(GSE3141_scpas_result$scRNA_data)[
  which(GSE3141_scpas_result$scRNA_data$scPAS == "Positive")
]
# GSE3141_scpas_neg <- colnames(GSE3141_scpas_result$scRNA_data)[
#     which(GSE3141_scpas_result$scRNA_data$scPAS == "Negative")
# ]

GSE3141_scab_result = qs::qread(
  "luad_GSE3141_scab.qs",
  nthreads = 4
)

GSE3141_scab_pos <- colnames(GSE3141_scab_result$scRNA_data)[
  which(GSE3141_scab_result$scRNA_data$scAB == "Positive")
]

# ! scPP is not suitable for this dataset
# GSE3141_scpp_result = qs::qread(
#     "luad_GSE3141_scpp.qs",
#     nthreads = 4
# )
# GSE3141_scpp_neg <- colnames(GSE3141_scpp_result$scRNA_data)[
#     which(GSE3141_scpp_result$scRNA_data$scPP == "Negative")
# ]

# ? ---- Sample info ----

all_cells = colnames(GSE3141_scissor_result$scRNA_data)

# > all_cells |>  length()
# [1] 32178

cnv_tumor_cells = colnames(GSE3141_scissor_result$scRNA_data)[which(
  GSE3141_scissor_result$scRNA_data$cnv_status == "tumor"
)]


# ? ----- TNBC tumor cell Venn diagram ----

pos_venn = list(
  scissor = GSE3141_scissor_pos,
  scpas = GSE3141_scpas_pos,
  scab = GSE3141_scab_pos,
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
    # "#9c8200ff",
    "blue",
    "purple",
    "cyan"
  )
) +
  ggplot2::scale_fill_gradient(low = "white", high = "#ffb6b6ff") +
  ggplot2::ggtitle("GSE3141 LUAD cells Venn diagram")

ggplot2::ggsave("GSE3141_LUAD_Venn.png", venn_plot, width = 10, height = 9)

# ? ----

seurat_merged = SigBridgeR::MergeResult(
  GSE3141_scab_result,
  GSE3141_scpas_result,
  GSE3141_scissor_result
)

qs::qsave(seurat_merged, "GSE3141_luad_merged_seurat.qs", nthreads = 4)
