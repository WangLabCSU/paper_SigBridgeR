data_path <- "/home/data/sigbridger/benchmark_data/brca"

GSE42568_bulk <- qs::qread(file.path(data_path, "brca_bulkdata_GSE42568.qs"))

genes <- rownames(GSE42568_bulk)
anyNA(genes)

unknown <- paste0("Unknown_", which(is.na(genes)))
genes[is.na(genes)] <- unknown
anyNA(genes)
