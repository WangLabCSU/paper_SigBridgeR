setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


set.seed(123)

library(dplyr)
library(data.table)
library(Seurat)

# ! tnbc-sc-GSE161529
# ! 阳性对照组
# ! binary phenotype

data_path = '/home/data/sigbridger/benchmark_binary/brca/TNBC'
seurat <- qs::qread(
  file.path(data_path, 'GSE42568_tnbc_merged_seurat.qs'),
  nthreads = 4L
)

screen_method <- c('scissor', 'scAB', 'scPP', 'scPAS')


for (method in screen_method) {
  if (!method %in% colnames(seurat[[]])) {
    warning('! ', method, ' not in seurat, skip it !')
    next
  }
  if (length(table(seurat[[method]])) < 2) {
    warning('! ', method, ' only have one class, skip it !')
    next
  }

  markers <- Seurat::FindAllMarkers(
    seurat,
    group.by = method,
    logfc.threshold = log(1.5), # 在两个不同组之间平均表达的差异倍数低于1.5的基因
    min.diff.pct = 0.1, # 过滤掉那些在两个不同组之间能检测到的细胞比例低于0.1的基因
    min.pct = 0.25 # 参数过滤掉那些在25%以下细胞中检测到的基因
  )

  markers_significant <- dplyr::filter(
    markers,
    p_val_adj < 0.05
  )

  only_top20 <- dplyr::group_by(markers_significant, cluster) %>%
    dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC))) %>%
    dplyr::slice_head(n = 20)

  data.table::fwrite(
    x = markers,
    file = paste0('binary_', method, '_sig_markers.csv'),
    row.names = TRUE
  )
  if (nrow(only_top20) < length(table(seurat[[method]])) * 20) {
    warning(
      '! ',
      method,
      ' one or more groups have less than 20 markers, maybe the threshold is too high, here skip it !'
    )
    next
  }

  data.table::fwrite(
    x = only_top20,
    file = paste0('binary_top20_', method, '_sig_markers.csv'),
    row.names = TRUE
  )
}
