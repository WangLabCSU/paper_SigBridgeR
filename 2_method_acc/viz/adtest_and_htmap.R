setwd(file.path(usethis::proj_path(), "2_method_acc/viz"))

library(dplyr)
library(data.table)
source("pairwise_ad_test.R")
source("ad_p_htmap.R")


combined2 <- data.table::fread("combined.csv")

# stats <- combined2[,
#   list(
#     N = .N,
#     mean_f1 = mean(F1, na.rm = TRUE),
#     sd_f1 = stats::sd(F1, na.rm = TRUE),
#     mean_acc = mean(Accuracy, na.rm = TRUE),
#     sd_acc = stats::sd(Accuracy, na.rm = TRUE),
#     median_f1 = stats::median(F1, na.rm = TRUE),
#     median_acc = stats::median(Accuracy, na.rm = TRUE),
#     min_f1 = min(F1, na.rm = TRUE),
#     max_f1 = max(F1, na.rm = TRUE),
#     min_acc = min(Accuracy, na.rm = TRUE),
#     max_acc = max(Accuracy, na.rm = TRUE)
#   ),
#   by = method_name
# ]

# stats2 <- combined2[,
#   list(
#     N = .N,
#     mean_f1 = mean(F1, na.rm = TRUE),
#     sd_f1 = stats::sd(F1, na.rm = TRUE),
#     mean_acc = mean(Accuracy, na.rm = TRUE),
#     sd_acc = stats::sd(Accuracy, na.rm = TRUE),
#     median_f1 = stats::median(F1, na.rm = TRUE),
#     median_acc = stats::median(Accuracy, na.rm = TRUE),
#     min_f1 = min(F1, na.rm = TRUE),
#     max_f1 = max(F1, na.rm = TRUE),
#     min_acc = min(Accuracy, na.rm = TRUE),
#     max_acc = max(Accuracy, na.rm = TRUE)
#   ),
#   by = list(method_name, tissue, screen_method)
# ]

ad_test_f1 <- pairwise_dist_test(
  dt = combined2,
  group_col = "method_name",
  value_col = "F1"
)
ad_test_acc <- pairwise_dist_test(
  dt = combined2,
  group_col = "method_name",
  value_col = "Accuracy"
)

data.table::fwrite(ad_test_f1, "ad_test_f1.csv")
data.table::fwrite(ad_test_acc, "ad_test_acc.csv")

ad_test_f1 <- data.table::fread("ad_test_f1.csv")
ad_test_acc <- data.table::fread("ad_test_acc.csv")

convert_to_p_mat <- function(wide) {
  wide <- tibble::as_tibble(wide)
  row_names <- unlist(wide$group1)
  total <- dplyr::union(colnames(wide[, -1]), row_names)

  missing <- dplyr::setdiff(total, row_names)
  missing_tbl <- data.frame(rep(NA, length(total)))
  missing_tbl <- t(missing_tbl)
  colnames(missing_tbl) <- total
  rownames(missing_tbl) <- missing
  res <- dplyr::bind_rows(
    tibble::column_to_rownames(wide, "group1"),
    as.data.frame(missing_tbl)
  )

  res <- res[sort(rownames(res)), sort(colnames(res))]
  res <- as.matrix(res)
  diag(res) <- 0L
  res <- Matrix::forceSymmetric(res, uplo = "U")
  as.matrix(res)
}

# * F1 score
wide_f1 <- ad_test_f1[, .(group1, group2, `AD.aympt. P-value`)]
wide_f1 <- data.table::dcast(
  wide_f1,
  formula = group1 ~ group2,
  value.var = "AD.aympt. P-value"
) %>%
  convert_to_p_mat()

# * Accuracy
wide_acc <- ad_test_acc[, .(group1, group2, `AD.aympt. P-value`)]
wide_acc <- data.table::dcast(
  wide_acc,
  formula = group1 ~ group2,
  value.var = "AD.aympt. P-value"
) %>%
  convert_to_p_mat()


# find_nonzero_min_val <- function(mat) {
#   mat[mat != 0 & !is.na(mat)] %>% min()
# }
# nearest_50 <- function(x) {
#   round(abs(x) / 50) * 50
# }

# tol_f1 <- find_nonzero_min_val(wide_f1) %>% log10() %>% nearest_50()
# tol_acc <- find_nonzero_min_val(wide_acc) %>% log10() %>% nearest_50()
# tol <- max(tol_f1, tol_acc)

tol <- 10L

ad_p_htmap(
  corr = wide_f1,
  p_mat = wide_f1,
  filename = "ad_p_htmap_f1.png",
  width = 3500,
  height = 3000,
  tolerance = 10^-tol
)


ad_p_htmap(
  corr = wide_acc,
  p_mat = wide_acc,
  filename = "ad_p_htmap_acc.png",
  width = 3500,
  height = 3000,
  tolerance = 10^-tol
)

cli::cli_h1("Done")

# ------------------------------------------------------------------------------------------
# legend
df <- data.frame(
  x = stats::rnorm(100),
  y = stats::rnorm(100),
  pval = stats::runif(100, 1e-8, 1)
)
df$neglogp <- -log10(df$pval)

# 颜色函数
col_fun <- grDevices::colorRampPalette(c(
  "#FEE5D9",
  "#FCAE91",
  "#FB6A4A",
  "#DE2D26",
  "#A50F15"
))

legend <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = neglogp)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_color_gradientn(
    colours = col_fun(100),
    name = "-log10(P value)",
    limits = c(0, tol)
  ) +
  ggplot2::theme(
    legend.title = ggplot2::element_text(face = "bold", size = 12)
  )

ggplot2::ggsave("adtest_and_htmap_legend.png", plot = legend)
