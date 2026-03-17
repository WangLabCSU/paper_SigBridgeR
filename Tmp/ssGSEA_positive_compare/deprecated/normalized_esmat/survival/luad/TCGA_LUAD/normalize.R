setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)

data_path <- "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA_positive_compare/esmat/survival/luad/TCGA_LUAD"
methods <- c("scissor", "scAB", "scPAS")

# ? read raw esmat
purrr::walk(
  methods,
  function(method_i) {
    esmat_file_path <- file.path(
      data_path,
      paste0(method_i, "_ssGSEA_score.qs")
    )

    assign(
      paste0("esmat_", method_i),
      qs::qread(esmat_file_path),
      envir = .GlobalEnv
    )
  },
  .progress = TRUE
)

# ? z-score normalize
purrr::walk(
  methods,
  function(method_i) {
    esmat_i <- get(paste0("esmat_", method_i))

    esmat_i <- esmat_i %>%
      group_by(cluster) %>%
      mutate(z_ssgsea_score = as.numeric(scale(ssgsea_score))) %>%
      ungroup()

    assign(
      paste0("esmat_", method_i),
      esmat_i,
      envir = .GlobalEnv
    )
  },
  .progress = TRUE
)

# ? save
purrr::walk(
  methods,
  function(method_i) {
    esmat_i <- get(paste0("esmat_", method_i))

    qs::qsave(
      esmat_i,
      file = paste0(method_i, "_ssGSEA_score_z.qs")
    )
  },
  .progress = TRUE
)
