data_path <- "/home/data/sigbridger/benchmark_data/brca"

GSE42568_bulk <- qs::qread(file.path(data_path, "brca_bulkdata_GSE42568.qs")) # transformed

GSE162228_bulk <- qs::qread(file.path(data_path, "brca_bulkdata_GSE162228.qs")) # transformed

TCGA_BRCA_bulk <- qs::qread(file.path(data_path, "brca_bulkdata_TCGA.qs")) # counts
