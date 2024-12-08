#' @include conda_S4_definitions.R
#' @include conda_S4_methods.R
NULL


################################################################################
##########################                            ##########################
#######                        class constructors                        #######
###                                                                          ###
#                                                                              #

#' Handy \code{CondaEnv} instance constructor
#'
#' @description
#' Wrapper to create a \linkS4class{CondaEnv} object meant to store
#' information about a single conda environment for a specific method. In most
#' cases, it is best that the user favors \code{CondaEnvManager}-related
#' functions.
#'
#' @param method the method name. One of "combat", "harmony", "mnn", "bbknn",
#' "scvi", "scanvi", "scanorama"
#' @param conda.bin the path to the conda binary. If empty (\code{""}),
#' \code{NULL} or \code{"conda"}, \code{"auto"} is passed to
#' \link[reticulate]{conda_binary}() to find the path to the conda binary in
#' \code{PATH}
#' @param conda.env.name name of the conda environment
#' @param conda.env.path path to the conda environment
#'
#' @return a \linkS4class{CondaEnv} object
#' @export
#' @seealso \linkS4class{CondaEnv}
CondaEnv <- function(method = known.methods, conda.bin = NULL,
                     conda.env.name = NULL, conda.env.path = NULL) {
  method <- tolower(method)
  method <- match.arg(method)
  conda.bin <- conda.bin %||% ''
  if (conda.bin == 'conda') conda.bin <- "auto"
  conda.env.name <- conda.env.name %||% ''
  conda.env.path <- conda.env.path %||% ''

  new("CondaEnv", method = method, needs.conda = needs.conda[[method]],
      conda.bin = new("CondaEnvSlot", value = conda.bin, valid = FALSE),
      conda.env.name = new("CondaEnvSlot", value = conda.env.name, valid = FALSE),
      conda.env.path = new("CondaEnvSlot", value = conda.env.path, valid = FALSE),
      is.valid = FALSE)
}

#' Handy \code{CondaEnvManger} instance constructor
#'
#' @description
#' Wrapper to create a \linkS4class{CondaEnvManager} object meant to store
#' information about conda environments for all known methods out of a data.frame
#' object (typically stored in the package's cache). In most cases, it is best
#' that the user favors the \link[=UpdateEnvCache]{UpdateEnvCache()} function
#' which enables to create/overwrite conda environments and to update the
#' package's internal cache.
#'
#' @param cache a correctly formatted (see below) data.frame or a path to such a
#' data.frame stored inside a RDS file. If cache is missing, return an blank
#' object (i.e. in an uninitialized state).
#'
#' @return a \linkS4class{CondaEnvManager} object
#'
#' @importFrom dplyr "%>%" group_by group_map
#' @export
#' @seealso \linkS4class{CondaEnvManager}
#' @details The \code{cache} data.frame's expected column names are:
#' \itemize{
#'   \item method
#'   \item needs.conda
#'   \item conda.bin.value
#'   \item conda.bin.valid
#'   \item conda.env.name.value
#'   \item conda.env.name.valid
#'   \item conda.env.path.value
#'   \item conda.env.path.valid
#'   \item is.valid
#' }
#'
CondaManager <- function(cache) {
  if (missing(cache)) {
    return(.blankCondaManager())
  } else {
    expected.fields <- c("method", "needs.conda", "conda.bin.value",
                         "conda.bin.valid",
                         "conda.env.name.value", "conda.env.name.valid",
                         "conda.env.path.value", "conda.env.path.valid",
                         "is.valid")
    if (is.character(cache)) {
      if (! file.exists(cache)) {
        abort(message = paste("provided path", sQuote(cache), "is not valid"))
      }
      cache <- readRDS(cache)
    }
    if (! all(expected.fields %in% colnames(cache)) || length(cache$method) != length(unique(cache$method))) {
      abort(message = "provided cache is not properly formatted (missing columns or non-unique methods)")
    }
    rownames(cache) <- cache$method
    missing.methods <- setdiff(known.methods, cache$method)
    unknown.methods <- setdiff(cache$method, known.methods)
    # if (length(missing.methods) == length(known.methods)) {
    #   warning("None of the cached method is known, initializing blank object",
    #           call. = TRUE, immediate. = TRUE)
    #   return(.blankCondaManager())
    # }
    if (length(unknown.methods)) {
      warning(paste(length(unknown.methods), " of the cached methods is/are unknown, skipping"),
              call. = TRUE, immediate. = TRUE)
    }
    cache.envs <- cache %>% group_by(method) %>%
      group_map(.f = .CondaEnvFromDf, .keep = TRUE)
    cache.envs <- append(cache.envs, lapply(missing.methods, CondaEnv))
    cache.envs <- setNames(cache.envs, sapply(cache.envs, function(obj) obj@method))
    conman <- suppressWarnings(
      new("CondaEnvManager",
          combat = checkCondaEnv(cache.envs$combat),
          harmony = checkCondaEnv(cache.envs$harmony),
          mnn = checkCondaEnv(cache.envs$mnn),
          bbknn = checkCondaEnv(cache.envs$bbknn),
          scvi = checkCondaEnv(cache.envs$scvi),
          scanvi = checkCondaEnv(cache.envs$scanvi),
          scanorama = checkCondaEnv(cache.envs$scanorama),
          trvae = checkCondaEnv(cache.envs$trvae))
    )
    if (length(missing.methods)) {
      warning(paste(length(missing.methods), " method(s) missing, adding"),
              call. = TRUE, immediate. = TRUE)
      saveToCache(conman)
    }
    return(conman)
  }
}

.blankCondaManager <- function() {
  suppressWarnings(
    new("CondaEnvManager",
        combat = checkCondaEnv(CondaEnv("combat")),
        harmony = checkCondaEnv(CondaEnv("harmony")),
        mnn = checkCondaEnv(CondaEnv("mnn")),
        bbknn = unsetSlot(CondaEnv("bbknn"), valid = FALSE),
        scvi = unsetSlot(CondaEnv("scvi"), valid = FALSE),
        scanvi = unsetSlot(CondaEnv("scanvi"), valid = FALSE),
        scanorama = unsetSlot(CondaEnv("scanorama"), valid = FALSE),
        trvae = unsetSlot(CondaEnv("trvae"), valid = FALSE))
  )
}

#' Handy \code{CondaEnvManger} instance modifier
#'
#' @description
#' Wrapper to update the \linkS4class{CondaEnvManager} object stored in cache
#' and used as configuration for conda environments for all known methods
#' (typically stored in the package's cache). A single method can be
#' updated at the time. In most cases, this is the best way to modify the
#' package's conda environments registry.
#'
#' @inheritParams CondaEnv
#' @param method the name of the method to update. One of "combat", "harmony",
#' "mnn", "bbknn", "scvi", "scanvi", "scanorama", "trvae"
#' @param conda.env Either the name of a conda environment or the path to such
#' an environment. Must be reachable by provided \code{conda.bin}. \code{NULL}
#' enable to use the default environment names. (see \strong{Details} section)
#' @param conda.env.is.path Whether the \code{conda.env} is a path (\code{TRUE}),
#' or a name (\code{FALSE}). The default ("\code{auto}") guesses. \strong{Be
#' careful not to make a mistake if you switch to a non-default value}
#' @param separate.scvi.envs By default, SCVI and SCANVI share the same conda
#' environment, since they rely on the same python package. If you wish to have
#' a distinct environment for each, set it to \code{TRUE}. Ignored if
#' \code{method} is not SCVI nor SCANVI.
#' @param overwrite.env Turn it to \code{TRUE} to enable overwriting an existing
#' environment (same name or same path). When the provided \code{conda.env}
#' already exists, the default behaviour is to update the package's registry
#' with the existing environment as it is
#' @param dry.run When \code{TRUE}, the package's current cache is not updated.
#' But the new conda environment (if any) will be created. \code{FALSE} by default.
#'
#' @return a \linkS4class{CondaEnvManager} object
#'
#' @export
#' @seealso \linkS4class{CondaEnvManager} \link[=CondaManager]{CondaManager()}
#' @details The conda environments default names are:
#' \itemize{
#'   \item \strong{bbknn}: \code{SeuratIntegrate_scvi}
#'   \item \strong{SCVI}: \code{SeuratIntegrate_scvi-tools}
#'   \item \strong{SCANVI}: \code{SeuratIntegrate_scvi-tools}
#'   \item \strong{scanorama}: \code{SeuratIntegrate_scanorama}
#'   \item \strong{trVAE}: \code{SeuratIntegrate_trvae}
#' }
#'
UpdateEnvCache <- function(method = known.methods, conda.bin = "auto",
                           conda.env = NULL, conda.env.is.path = "auto",
                           separate.scvi.envs = FALSE, overwrite.env = FALSE,
                           dry.run = FALSE) {
  conman <- conda_status$current
  method <- tolower(method)
  method <- match.arg(method)

  if (! needs.conda[[method]]) {
    message(sprintf("%s does not require a conda environment. Skipping", method))
    return(invisible(0))
  }

  conda.bin <- conda.bin %||% 'auto'
  if (conda.bin %in% c('', 'conda')) {
    conda.bin <- 'auto'
  }
  conda.env %||% {conda.env <- .conda_default_envname[[method]]
                  conda.env.is.path <- FALSE}
  conda.env.is.path <- conda.env.is.path %||% "auto"
  if (conda.env.is.path == "auto") {
    conda.env.is.path <- grepl("^~|/|\\\\", conda.env) |
      (dir.exists(dirname(conda.env)) & dirname(conda.env) != ".")
  }
  conda.env.name <- conda.env
  conda.env.path <- NULL
  if (conda.env.is.path) {
    conda.env.path <- conda.env <- normalizePath(conda.env)
    conda.env.name <- basename(conda.env.path)
  }
  separate.scvi.envs <- separate.scvi.envs %||% FALSE

  conda.known.envs <- conda_list(conda = conda.bin)
  conda.known.envs$path <- sub("/bin/python[^/]*$", "", conda.known.envs$python)

  create.env <- TRUE
  if (condaenv_exists(envname = conda.env, conda = conda.bin)) {
    if (overwrite.env) {
      message("environment ", sQuote(conda.env), " already exists. ",
              "Removing and setting fresh new environment")
      conda_remove(conda.env, packages = NULL, conda = conda.bin)
    } else {
      message("environment ", sQuote(conda.env), " already exists. ",
              "Updating known conda environments and leaving")
      create.env <- FALSE
    }
  }
  if (create.env) {
    conda.env_ <- conda.env
    os <- ifelse("Darwin" %in% Sys.info()[['sysname']], "osx", "default")
    conda.req <- .conda_requirements[[os]][[method]]
    conda.env <- conda_create(envname = conda.env_,
                              packages = conda.req$packages,
                              forge = FALSE,
                              channel = conda.req$channels,
                              conda = conda.bin)
    conda.known.envs <- conda_list(conda = conda.bin)
    conda.known.envs$path <- sub("/bin/python[^/]*$", "", conda.known.envs$python)

    conda.req$pip %iff% {
      conda_install(packages = conda.req$pip,
                    envname = conda.env_,
                    pip = TRUE, conda = conda.bin)
    }
  }
  if (conda.env.is.path) {
    conda.env.name <- conda.known.envs$name[which.min(adist(conda.env, conda.known.envs$path))]
  } else {
    conda.env.path <- conda.known.envs$path[conda.known.envs$name == conda.env.name][1]
  }
  conman <- updateEnv(conman,
                      method = method,
                      conda.bin = conda.bin,
                      conda.env.name = conda.env.name,
                      conda.env.path = conda.env.path)
  if (method %in% c("scvi", "scanvi")) {
    other.scvi <- setdiff(c("scvi", "scanvi"), method)
    if (! separate.scvi.envs) {
      slot(object = conman, name = other.scvi) <- slot(object = conman, name = method)
      slot(object = conman, name = other.scvi)@method <- other.scvi
    }# else {
    #   slot(object = conman, name = other.scvi) <- CondaEnv(other.scvi)
    # }
  }

  if (!dry.run) {
    conda_status$current <- conman
    saveToCache(conman)
    return(invisible(conman))
  } else {
    return(conman)
  }
}

.CondaEnvFromDf <- function(cache, ...) {
  new("CondaEnv", method = cache$method, needs.conda = cache$needs.conda,
      conda.bin = new("CondaEnvSlot",
                      value = cache$conda.bin.value,
                      valid = cache$conda.bin.valid),
      conda.env.name = new("CondaEnvSlot",
                           value = cache$conda.env.name.value,
                           valid = cache$conda.env.name.valid),
      conda.env.path = new("CondaEnvSlot",
                           value = cache$conda.env.path.value,
                           valid = cache$conda.env.path.valid),
      is.valid = cache$is.valid)
}

CondaEnvFromDf <- function(cache) {
  expected.fields <- c(
    "method", "needs.conda", "conda.bin.value", "conda.bin.valid",
    "conda.env.name.value", "conda.env.name.valid", "conda.env.path.value",
    "conda.env.path.valid", "is.valid")
  if (! all(expected.fields %in% names(cache))) {
    abort("Bad cache formatting")
  }
  .CondaEnvFromDf(cache)
}

.blankCondaManager <- function() {
  new("CondaEnvManager",
      combat = checkCondaEnv(CondaEnv("combat")),
      harmony = checkCondaEnv(CondaEnv("harmony")),
      mnn = checkCondaEnv(CondaEnv("mnn")),
      bbknn = unsetSlot(CondaEnv("bbknn"), valid = FALSE),
      scvi = unsetSlot(CondaEnv("scvi"), valid = FALSE),
      scanvi = unsetSlot(CondaEnv("scanvi"), valid = FALSE),
      scanorama = unsetSlot(CondaEnv("scanorama"), valid = FALSE),
      trvae = unsetSlot(CondaEnv("trvae"), valid = FALSE))
}

#' Get path to package config cache
#'
#' @description
#' Get path to package config cache file or directory
#'
#' @param include.file Set to \code{FALSE} to retrieve the directory containing
#' the cache file. Return the path to the file by default.
#'
#' @return A single character, the path to the cache file or directory
#'
#' @importFrom tools R_user_dir
#' @export
getCachePath <- function(include.file = TRUE) {
  include.file <- ifelse(is.logical(include.file %||% FALSE), include.file[[1]], TRUE)
  file.path(R_user_dir(.packageName, which = 'config'),
            c('', 'conda_envs.RDS')[include.file + 1])
}

#' Get current cache of conda environments
#'
#' @description
#' Get current cache of conda environments from package environment or cache on disk
#'
#' @param from Either "env" or "cache". Where to load cache from.
#'
#' @return a \linkS4class{CondaEnvManager} object
#'
#' @export
#'
getCache <- function(from = c("env", "cache")) {
  from <- tolower(from)
  from <- match.arg(from)
  if (from == "env") {
    return(conda_status$current)
  }
  return(CondaManager(getCachePath(include.file = TRUE)))
}

#' Reload cache from disk
#'
#' @description
#' Reload cache from disk to update current list of conda environments in
#' package environment
#'
#' @return a \linkS4class{CondaEnvManager} object (invisibly)
#'
#' @export
reloadCache <- function() {
  conda_status$current <- getCache(from = "cache")
  return(invisible(conda_status$current))
}

#' Reset a conda environment
#'
#' @description
#' Reset (unset) a conda environment linked to a provided method. Update cache
#' on disk and in package environment
#' @inheritParams UpdateEnvCache
#'
#' @return a \linkS4class{CondaEnvManager} object (invisibly)
#'
#' @export
resetCache <- function(method = known.methods) {
  method <- tolower(method)
  method <- match.arg(method)

  conda_status$current[[method]] <- CondaEnv(method = method,
                                             conda.bin = NULL,
                                             conda.env.name = NULL,
                                             conda.env.path = NULL)
  saveToCache(conda_status$current)
  return(invisible(conda_status$current))
}
