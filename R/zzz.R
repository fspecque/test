#' @export
conda_status <- new.env(parent = emptyenv())

pb <- new.env(parent = parent.frame())

.onLoad <- function(libname = find.package(.packageName), pkgname = .packageName) {
  cache.path <-  getCachePath(include.file = TRUE)
  pb$onLoad <- TRUE
  pb$onAttach <- FALSE
  if (!file.exists(cache.path)) {
    cache.path <- NULL
  }
  conda_status$current <- CondaManager(cache.path)
  pb$onLoad <- FALSE
  pb$print.lisi.msg <- TRUE
}

.onAttach <- function(libname = find.package(.packageName), pkgname = .packageName) {
  cache.path <- getCachePath(include.file = FALSE)
  if (!dir.exists(cache.path)) {
    dir.create(cache.path, recursive = TRUE)
    packageStartupMessage("Cache directory created at ", cache.path)
  }

  cache.path <-  getCachePath(include.file = TRUE)
  invisible(utils::capture.output(
    pb$pb <- utils::txtProgressBar(0, sum(needs.conda), width = 30, style = 3)))
  pb$onLoad <- FALSE
  pb$onAttach <- TRUE
  w <- "found"
  if (!file.exists(cache.path)) {
    saveToCache(.blankCondaManager())
    w <- "created"
  }
  packageStartupMessage("Cache file ", w , " at ", cache.path, ". Loading...")
  conda_status$current <- CondaManager(cache.path)
  packageStartupMessage(utils::capture.output(close(pb$pb)))
  pb$onAttach <- FALSE
  pb$print.lisi.msg <- TRUE
}
