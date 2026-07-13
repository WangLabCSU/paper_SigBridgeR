setwd(dir = file.path(usethis::proj_path(), "5_GSEA"))
geneset_hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
geneset_hallmark <- data.table::as.data.table(geneset_hallmark)
data.table::fwrite(x = geneset_hallmark, file = "geneset_hallmark.csv")

geneset_desc <- split(geneset_hallmark$gs_description, geneset_hallmark$gs_name)
names(geneset_desc) <- stringr::str_remove(names(geneset_desc), "HALLMARK_")

geneset_desc <- lapply(geneset_desc, base::unique)
data.table::fwrite(x = geneset_desc, file = "geneset_hallmark_desc.csv")
