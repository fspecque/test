#' @export
conda_status <- new.env(parent = emptyenv())

pb <- new.env(parent = emptyenv())

.onLoad <- function(libname = find.package(.packageName), pkgname = .packageName) {
  cache.path <- getCachePath(include.file = FALSE)
  if (!dir.exists(cache.path)) {
    dir.create(cache.path)
    message("Cache directory created at ", cache.path)
  }

  cache.path <-  getCachePath(include.file = TRUE)
  w <- "found"
  if (!file.exists(cache.path)) {
    saveToCache(.blankCondaManager())
    w <- "created"
  }
  message("Cache file ", w , " at ", cache.path, ". Loading...")
  invisible(capture.output(
    pb$pb <- utils::txtProgressBar(0, sum(needs.conda), width = 30, style = 3)))
  pb$onLoad <- TRUE
  conda_status$current <- CondaManager(cache.path)
  pb$onLoad <- FALSE
}
