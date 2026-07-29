setwd(usethis::proj_path())

library(UCSCXenaShiny)
data(tcga_clinical)
data(tcga_clinical_fine)
data(tcga_surv)

source("0_bench_data_preprocess/row_names_duplifinder.R")

out_dir <- "/home/data/sigbridger/benchmark_data/"

table(tcga_clinical$type)

luad <- dplyr::filter(tcga_clinical, type == "LUAD") %>%
  dplyr::mutate(tissue_type = stringr::str_remove(sample, ".*-")) %>%
  dplyr::filter(tissue_type %in% c("01", "11")) %>%
  dplyr::mutate(
    tumor = dplyr::case_when(tissue_type == "01" ~ 1, tissue_type == "11" ~ 0)
  )

FindDuplicates(luad, luad$sample)

luad <- dplyr::distinct(luad, sample, .keep_all = TRUE)

brca <- dplyr::filter(tcga_clinical, type == "BRCA") %>%
  dplyr::mutate(tissue_type = stringr::str_remove(sample, ".*-")) %>%
  dplyr::filter(tissue_type %in% c("01", "11")) %>%
  dplyr::mutate(
    tumor = dplyr::case_when(tissue_type == "01" ~ 1, tissue_type == "11" ~ 0)
  )

FindDuplicates(brca, brca$sample)

rownames(brca) <- brca$sample
rownames(luad) <- luad$sample

qs::qsave(brca, file.path(out_dir, "brca/brca_pheno_TCGA.qs"))
qs::qsave(luad, file.path(out_dir, "lung/TCGA_LUAD_pheno.qs"))

# ------------------------------------------------------------------------------------------------------

luad_samples <- luad$sample
brca_samples <- brca$sample
luad_surv <- dplyr::filter(tcga_surv, sample %in% luad_samples) %>%
  dplyr::distinct(sample, .keep_all = TRUE) %>%
  tibble::column_to_rownames("sample") %>%
  dplyr::select(OS.time, OS) %>%
  dplyr::rename(time = 1, status = 2)
brca_surv <- dplyr::filter(tcga_surv, sample %in% brca_samples) %>%
  dplyr::distinct(sample, .keep_all = TRUE) %>%
  tibble::column_to_rownames("sample") %>%
  dplyr::select(OS.time, OS) %>%
  dplyr::rename(time = 1, status = 2)

qs::qsave(luad_surv, file.path(out_dir, "lung/TCGA_LUAD_surv_pheno.qs"))
qs::qsave(brca_surv, file.path(out_dir, "brca/brca_surv_TCGA.qs"))
