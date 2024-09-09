devtools::document()
pkgload::dev_help("AddLISIScore")

DESCRIPTION <- readLines("DESCRIPTION")
pkg.version <- gsub("\\s*", "", sub("^Version:", "", DESCRIPTION[grepl("^Version", DESCRIPTION)]))
archive.name <- paste0('SeuratIntegrate_', pkg.version,'.tar.gz')
cmd = paste0("R CMD build --no-build-vignettes ../Seurat-Integrate && mv ", archive.name,
             " ../ && R CMD INSTALL ../", archive.name)
system(cmd)
