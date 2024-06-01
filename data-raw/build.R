# roxygen2::roxygenise()
devtools::document()
pkgload::dev_help("ScoreRegressPC.CellCycle")
# pkgload::dev_help("ScoreDensityPC")

DESCRIPTION <- readLines("DESCRIPTION")
pkg.version <- gsub("\\s*", "", sub("^Version:", "", DESCRIPTION[grepl("^Version", DESCRIPTION)]))
archive.name <- paste0('SeuratIntegrate_', pkg.version,'.tar.gz')
cmd = paste0("R CMD build ../Seurat-Integrate && mv ", archive.name,
             " ../ && R CMD INSTALL ../", archive.name)
system(cmd)
