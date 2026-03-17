setwd('/home/yyx/R/Project/R_code/SigBridgeR')
devtools::document("..")

data_path = '/home/data/sigbridger'
source('/home/yyx/R/Project/R_code/SigBridgeR/R/90-Test.R')


scissor_result <- qs::qread(
  file.path(data_path, 'scissor_result.qs'),
  nthreads = 4L
)
scab_result <- qs::qread(
  file.path(data_path, 'scab_result.qs'),
  nthreads = 4L
)
scpas_result <- qs::qread(
  file.path(data_path, 'scpas_result.qs'),
  nthreads = 4L
)
scpp_result <- qs::qread(
  file.path(data_path, 'scpp_result.qs'),
  nthreads = 4L
)
degas_result <- qs::qread(
  file.path(data_path, 'degas_result.qs'),
  nthreads = 4L
)
lpsgl_result <- qs::qread(
  file.path(data_path, 'lpsgl_result.qs'),
  nthreads = 4L
)
# ? mistaken names
names(lpsgl_result)[1] <- "scRNA_data"


# ? binary phenotype
pipet_result <- qs::qread(
  file.path(data_path, 'pipet_result.qs'),
  nthreads = 4L
)

merged_without_pipet <- MergeResult(
  scissor_result,
  scab_result,
  scpas_result,
  scpp_result,
  degas_result,
  lpsgl_result
)

qs::qsave(
  merged_without_pipet,
  file.path(data_path, 'merged_without_pipet.qs'),
  nthreads = 4L
)


c(scissor_pos, scab_pos, scpas_pos, scpp_pos, degas_pos, lpsgl_pos) %<-%
  purrr::map(
    c("scissor", "scAB", "scPAS", "scPP", "DEGAS", "LP_SGL"),
    ~ colnames(merged_without_pipet)[
      which(merged_without_pipet[[.x]] == "Positive")
    ]
  )

all_cells <- colnames(merged_without_pipet)


ggplot2::ggsave(
  "vignettes/example_figures/venn.png",
  plot = venn,
  width = 10,
  height = 10
)


upset <- ScreenUpset(
  screened_seurat = merged_without_pipet,
  screen_type = c("scissor", "scPAS", "scAB", "scPP", "DEGAS", "LP_SGL"),
  n_intersections = 40
)

ggplot2::ggsave(
  "vignettes/example_figures/upset.png",
  plot = upset$plot,
  width = 10,
  height = 10
)

set.seed(123)
# * fictional sample column
merged_without_pipet$Sample <- sample(
  paste0("Sample", 1:10),
  ncol(merged_without_pipet),
  replace = TRUE
)


fraction_list = ScreenFractionPlot(
  screened_seurat = merged_without_pipet,
  group_by = "Sample",
  screen_type = c("scissor", "scPP", "scAB", "scPAS", "DEGAS", "LP_SGL"),
  show_null = FALSE,
  plot_color = NULL,
  show_plot = TRUE
)

ggplot2::ggsave(
  "vignettes/example_figures/fraction.png",
  plot = fraction_list$combined_plot,
  width = 10,
  height = 10
)


my_palette <- randomcoloR::distinctColorPalette(
  k = length(unique(merged_without_pipet$Sample)),
  runTsne = TRUE
)

sample_umap = Seurat::DimPlot(
  merged_without_pipet,
  group.by = "Sample",
  pt.size = 1.2,
  alpha = 0.8,
  reduction = "umap",
  cols = my_palette
) +
  ggplot2::ggtitle("Sample")


c(
  scissor_umap,
  scab_umap,
  scpas_umap,
  scpp_umap,
  degas_umap,
  lpsgl_umap
) %<-%
  purrr::map(
    c("scissor", "scAB", "scPAS", "scPP", "DEGAS", "LP_SGL"), # make sure these column names exist
    ~ Seurat::DimPlot(
      merged_without_pipet,
      group.by = .x,
      pt.size = 1.2,
      alpha = 0.8,
      reduction = "umap",
      cols = c(
        "Neutral" = "#CECECE",
        "Other" = "#CECECE",
        "Positive" = "#ff3333",
        "Negative" = "#386c9b"
      )
    ) +
      ggplot2::ggtitle(.x)
  )


ggplot2::ggsave(
  "vignettes/example_figures/umaps.png",
  plot = umaps,
  width = 10,
  height = 10
)


set.seed(123)

gene_ids <- c(
  paste0("ENSG", sprintf("%011d", 1:3900)),
  paste0("GENE_", LETTERS[1:100])
)
library(MASS)
mu <- c(
  runif(500, 100, 2000), # high
  runif(1500, 20, 200), # medium
  rexp(2000, rate = 1 / 2) # low (gamma-like, many near-zero)
)
lib_size_factor <- rlnorm(506, meanlog = log(3e7), sdlog = 0.3) / 3e7 # ~20–60M reads
counts_mat <- matrix(0L, nrow = 4000, ncol = 506)

for (i in seq_len(4000)) {
  mu_i <- mu[i] * lib_size_factor
  # dispersion decreases with mean: phi ≈ 0.5 / sqrt(mu_i + 1)
  phi_i <- pmax(0.05, 0.5 / sqrt(mu_i + 1))
  size_i <- 1 / phi_i
  counts_mat[i, ] <- rnbinom(506, mu = mu_i, size = size_i)
}
counts_mat <- floor(counts_mat)
counts_mat[counts_mat < 0] <- 0L

colnames(counts_mat) <- sample_info$sample
rownames(counts_mat) <- gene_ids
