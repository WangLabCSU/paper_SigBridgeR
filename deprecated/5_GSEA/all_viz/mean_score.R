# ! all tissue

setwd(file.path(usethis::proj_path(), "5_GSEA/all_viz"))

library(dplyr)
library(data.table)

data_dir <- "/home/data/sigbridger/GSEA"
data_files <- list.files(
  data_dir,
  pattern = "binary.*\\.csv|survival.*\\.csv",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE,
  include.dirs = FALSE
)
names(data_files) <- basename(data_files) %>%
  tools::file_path_sans_ext() %>%
  gsub("_merge.*", "", .)

# * "Negative""Neutral" -> "Other"
mirai::daemons(4L)
loaded_data <- mirai::mirai_map(data_files, \(data_path) {
  library(data.table)
  method_cols <- c(
    "scissor",
    "scAB",
    "scPAS",
    "scPP",
    "DEGAS",
    "PIPET",
    "SCIPAC",
    "LP_SGL"
  )
  loaded <- data.table::fread(data_path)
  existing_cols <- colnames(loaded)[colnames(loaded) %in% method_cols]

  loaded[,
    (existing_cols) := lapply(.SD, \(COL) {
      data.table::fifelse(COL != "Positive", "Other", "Positive")
    }),
    .SDcols = existing_cols
  ]
})
loaded_data <- loaded_data[mirai::.progress]
mirai::daemons(0L)

# * save intermediate data
if (!file.exists(file.path(data_dir, "all_gsea_score.qs"))) {
  qs::qsave(
    loaded_data,
    file.path(data_dir, "all_gsea_score.qs"),
    nthreads = 8L
  )
  cli::cli_alert_success("Saved {.path all_gsea_score.qs}")
} else {
  cli::cli_warn("all_gsea_score.qs already exists!")
}

# * stats by group
mean_score_gsea <- lapply(loaded_data, \(dt) {
  col_names <- colnames(dt)
  existing_cols <- col_names[col_names %in% method_cols]

  # gsea score
  value_cols <- col_names[grepl("^A|^Gene", col_names)]
  # delete unused cols
  dt[, (setdiff(colnames(dt), c(value_cols, existing_cols))) := NULL]

  # 一步到位：melt → 分组计数 + 求均值
  data.table::melt(
    dt,
    id.vars = value_cols,
    measure.vars = existing_cols,
    variable.name = "method",
    value.name = "status"
  )[,
    c(
      .(count = .N),
      lapply(.SD, mean, na.rm = TRUE)
    ),
    by = .(method, status),
    .SDcols = value_cols
  ]
})

dir.create(
  file.path(data_dir, "mean_score_gsea"),
  recursive = TRUE,
  showWarnings = FALSE
)
purrr::iwalk(
  mean_score_gsea,
  \(dt, dt_name) {
    data.table::fwrite(
      dt,
      file.path(
        data_dir,
        "mean_score_gsea",
        paste0(dt_name, "_mean_hallmark_score.csv")
      )
    )
  },
  .progress = "Writing"
)
