library(dplyr)
library(data.table)
library(ggplot2)

# ! scissor does not need iteration, this script is deprecated.

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/cross_args_ite/GSE165897/GSE140082"
)

data_path = "/home/data/sigbridger/benchmark_data/ov"
devtools::document("~/R/Project/R_code/SigBridgeR")

ncore <- 16 # 留 1 核给系统
cat("使用核心数:", ncore, "\n")

# * load data
seurat = qs::qread(file.path(data_path, "hgsoc_GSE165897_seurat.qs"))

bulk = qs::qread(
  file.path(data_path, "ov_bulkdata_GSE140082.qs")
)
pheno = qs::qread(file.path(data_path, "ov_pheno_GSE140082.qs"))

surv_data = select(pheno, "OS (M)", "Death (1)", "Sample_ID") %>%
  filter(Sample_ID %in% colnames(bulk)) %>%
  rename("time" = 1, "status" = 2) %>%
  tibble::column_to_rownames("Sample_ID") %>%
  mutate_all(~ as.numeric(.))

# ! -------------------- scissor ---------------------

# * random search, 100 times
set.seed(123)


# * run scissor
if (.Platform$OS.type == "unix") {
  res_list <- parallel::mclapply(
    seq_len(nrow(alpha_samples)),
    function(i) {
      scissor_result = Screen(
        matched_bulk = bulk,
        sc_data = seurat,
        phenotype = pheno,
        label_type = glue::glue("OS (M)_survival_{i}"),
        phenotype_class = "survival",
        screen_method = "Scissor"
      )

      qs::qsave(
        scissor_result,
        file = glue::glue("scissor_results/scissor_result_{i}.qs")
      )

      pos_ratio = mean(scissor_result$scRNA_data$Scissor == "Positive")

      return(pos_ratio)
    },
    mc.cores = ncore
  )
}

# *visualize
gc()
res_list = unlist(res_list)
alpha_samples$pos_ratio = res_list

plot = ggplot(alpha_samples, aes(x = alpha, y = alpha2, color = pos_ratio)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "white", high = "red", limits = c(0, 1)) +
  theme_bw() +
  ggtitle(
    "Validation of the Screening Efficiency of scissor under Random Parameters"
  )

write.csv(alpha_samples, file = "scissor_random_search.csv")

ggsave(
  filename = "scissor_random_search.png",
  plot = plot,
  width = 8,
  height = 8,
  dpi = 300
)
