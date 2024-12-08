#' @include conda_S4_definitions.R
#' @importFrom reticulate conda_binary condaenv_exists conda_list conda_create conda_remove conda_install
NULL

################################################################################
##########################                            ##########################
#######                           set generics                           #######
###                                                                          ###
#                                                                              #

#' Check the validity of conda environment's components
#'
#' @description
#' Check the validity of conda environment's components. Best used on a
#' \linkS4class{CondaEnv} object. Internally call:
#' \itemize{
#'   \item \strong{checkCondaBin}: \link[reticulate]{conda_binary}
#'   \item \strong{checkCondaEnvPath}: \link[reticulate]{condaenv_exists}
#'   \item \strong{checkCondaEnvName}: \link[reticulate]{condaenv_exists}
#'   \item \strong{checkCondaEnv}: both
#' }
#'
#' @param x a \linkS4class{CondaEnv} object or a sub element of it of class
#' \code{CondaEnvSlot}. \code{checkCondaEnv()} only accepts a CondaEnv object
#' @param ... ignored when \code{x} is a \linkS4class{CondaEnv} object.
#' Otherwise, additional arguments to \pkg{reticulate} functions in
#' \strong{Description} such as \code{conda = ...}
#' @param verbose ignored
#'
#' @return an updated object of same class as input \code{x}
#' @export
#' @note \code{checkCondaEnv} is a wrapper that calls all checks.
#' @seealso \link[=isValid]{isValid()}, \linkS4class{CondaEnv}
#' @rdname checkCondaEnv
setGeneric("checkCondaBin", function(x, ..., verbose = getOption("verbose"))
  standardGeneric("checkCondaBin"))
#' @export
#' @rdname checkCondaEnv
setGeneric("checkCondaEnvPath", function(x, ..., verbose = getOption("verbose"))
  standardGeneric("checkCondaEnvPath"))
#' @export
#' @rdname checkCondaEnv
setGeneric("checkCondaEnvName", function(x, ..., verbose = getOption("verbose"))
  standardGeneric("checkCondaEnvName"))
#' @export
#' @rdname checkCondaEnv
setGeneric("checkCondaEnv", function(x, ..., verbose = getOption("verbose"))
  standardGeneric("checkCondaEnv"))

#' Check the validity of conda environment's components
#'
#' @description
#' Get the validity of conda environment's components. Best used on a
#' \linkS4class{CondaEnvManager} or a
#' \linkS4class{CondaEnv} object.
#'
#' @param x a \linkS4class{CondaEnvManager}, a
#' \linkS4class{CondaEnv} object or a sub element of it of class
#' \code{CondaEnvSlot}.
#' @param do.check whether to call \link[=checkCondaEnv]{check functions}. Does
#' not apply when \code{x} is a \code{CondaEnvSlot} object.
#' @param ... ignored when \code{x} is a \code{CondaEnvSlot} object. Otherwise,
#' can be used to pass \code{do.check} argument.
#'
#' @return \code{TRUE} or \code{FALSE} or a logical named vector when \code{x}
#' is a \linkS4class{CondaEnvManager} object (one element per
#' method)
#' @export
#' @seealso \linkS4class{CondaEnvManager}, \linkS4class{CondaEnv},
#' \link[=checkCondaEnv]{check conda functions}
setGeneric("isValid", function(x, ...)
  standardGeneric("isValid"))
setGeneric("unsetSlot", function(x, ...)
  standardGeneric("unsetSlot"))

#' Save a conda environment manager to cache
#'
#' @description
#' Save a conda environment manager to cache (
#' \code{\link[tools:R_user_dir]{R_user_dir("SeuratIntegrate", which = "config")}})
#' as a data.frame in a compressed RDS file.
#'
#' @param x a \linkS4class{CondaEnvManager} object
#' @param ... ignored
#' @param verbose ignored
#'
#' @export
#' @seealso \linkS4class{CondaEnvManager}
setGeneric("saveToCache", function(x, ..., verbose = getOption("verbose"))
  standardGeneric("saveToCache"))
setGeneric("updateEnv", function(x, ...)
  standardGeneric("updateEnv"))


################################################################################
##########################                            ##########################
#######                            set methods                           #######
###                                                                          ###
#                                                                              #

# Dispatch methods
setMethod("checkCondaBin", "CondaEnvSlot", function(x) {
  bin.path <- tryCatch(conda_binary(x@value),
                       error = function(e) NULL)
  x@value <- bin.path %||% x@value
  x@valid <- bin.path %iff% TRUE %||% FALSE
  x
})
setMethod("checkCondaBin", "CondaEnv", function(x) {
  if (x@needs.conda) {
    x@conda.bin <- checkCondaBin(x@conda.bin)
  } else {
    x <- unsetSlot(x = x, slot = "conda.bin", valid = TRUE)
  }
  x
})

setMethod("checkCondaEnvPath", "CondaEnvSlot", function(x, ...) {
  env.exists <- condaenv_exists(x@value, ...) & file.exists(x@value)
  x@valid <- env.exists
  if (env.exists) {
    x@value <- normalizePath(x@value)
  }
  x
})
setMethod("checkCondaEnvPath", "CondaEnv", function(x, ...) {
  if (x@needs.conda) {
    x@conda.env.path <- checkCondaEnvPath(x@conda.env.path,
                                          conda = x@conda.bin@value)
  } else {
    x <- unsetSlot(x = x, slot = "conda.env.path", valid = TRUE)
  }
  x
})

setMethod("checkCondaEnvName", "CondaEnvSlot", function(x, ...) {
  env.exists <- condaenv_exists(x@value, ...)
  x@valid <- env.exists
  x
})
setMethod("checkCondaEnvName", "CondaEnv", function(x, ...) {
  if (x@needs.conda) {
    x@conda.env.name <- checkCondaEnvName(x@conda.env.name,
                                          conda = x@conda.bin@value)
  } else {
    x <- unsetSlot(x = x, slot = "conda.env.name", valid = TRUE)
  }
  x
})

#' @importFrom utils setTxtProgressBar
setMethod("checkCondaEnv", "CondaEnv", function(x, ...) {
  if (pb$onLoad) {
    tryCatch({
      setTxtProgressBar(pb$pb, which(known.methods[needs.conda] == x@method))
    }, error = function(e) NULL)
  }
  x <- checkCondaBin(x)
  x <- checkCondaEnvPath(x)#, conda = x@conda.bin@value)
  x <- checkCondaEnvName(x)#, conda = x@conda.bin@value)
  if (! isValid(x@conda.env.path) & isValid(x@conda.env.name)) {
    conda.list.df <- conda_list(conda = x@conda.bin@value)
    x@conda.env.path@value <- basename(
      sub("/bin/python[^/]*$", "",
          conda.list.df$python[conda.list.df$name == x@conda.env.name@value]))
    x <- checkCondaEnvPath(x)#, conda = x@conda.bin@value)
    if (! isValid(x@conda.env.path)) {
      warning("Could not deduce path to conda environment base on name. Invalid",
              call. = T, immediate. = T)
    }
  } else if (isValid(x@conda.env.path) & ! isValid(x@conda.env.name)) {
    x@conda.env.name@value <- basename(x@conda.env.path@value)
    x <- checkCondaEnvName(x, conda = x@conda.bin)
    if (! isValid(x@conda.env.name)) {
      warning("Could not deduce name of conda environment base on path. Invalid",
              call. = T, immediate. = T)
    }
  }
  if (x@conda.env.name@value != basename(x@conda.env.path@value)) {
    x@conda.env.name@valid <- FALSE
    x@conda.env.path@valid <- FALSE
    warning(sprintf("provided conda environment name (%s) does not match provided path to conda environment (%s)",
                    sQuote(x@conda.env.name@value), sQuote(x@conda.env.path@value)),
            call. = T, immediate. = T)
  }
  are.valid <- sapply(c("conda.bin", "conda.env.name", "conda.env.path"),
                      FUN = slot, object = x, simplify = F)
  are.valid <- sapply(are.valid, FUN = isValid, simplify = "array")
  # x@is.set <- all(are.valid)
  x@is.valid <- all(are.valid)
  x
})
setMethod("checkCondaEnv", "CondaEnvManager", function(x, ...) {
  list.of.envs <- lapply(x, checkCondaEnv)
  invisible(lapply(names(list.of.envs), function(env) {
    slot(object = x, name = env) <<- list.of.envs[[env]]
  }))
  x
})

#' @export
#' @rdname isValid
setMethod("isValid", "CondaEnvSlot", function(x, ...) {
  x@valid
})
#' @export
#' @rdname isValid
setMethod("isValid", "CondaEnv", function(x, do.check = FALSE, ...) {
  if (do.check) {
    x <- checkCondaEnv(x)
  }
  x@is.valid
})
#' @export
#' @rdname isValid
setMethod("isValid", "CondaEnvManager", function(x, ...) {
  sapply(x, isValid, ...)
})

setMethod("unsetSlot", "CondaEnvSlot", function(x, valid = FALSE) {
  x@value <- ''
  x@valid <- valid %||% FALSE
  x
})

#' @importFrom methods slot slot<-
setMethod("unsetSlot", "CondaEnv", function(x, slot = c('conda.bin',
                                                        'conda.env.name',
                                                        'conda.env.path'),
                                            valid = FALSE) {
  slot <- tolower(slot)
  slot <- match.arg(slot, several.ok = TRUE)
  for (name in slot) {
    slot(object = x, name = name) <- unsetSlot(slot(object = x, name = name),
                                               valid = valid)
  }
  x
})

setMethod("saveToCache", "CondaEnvManager", function(x) {
  cache <- as.data.frame(x)
  saveRDS(cache, getCachePath(include.file = TRUE),
          ascii = FALSE, compress = "gzip")
})

setMethod("updateEnv", "CondaEnv", function(x,
                                            conda.bin = NULL,
                                            conda.env.name = NULL,
                                            conda.env.path = NULL) {
  new.env <- checkCondaEnv(CondaEnv(method = x@method, conda.bin = conda.bin,
                                    conda.env.name = conda.env.name,
                                    conda.env.path = conda.env.path))
  new.env
})
setMethod("updateEnv", "CondaEnvManager", function(x, method = known.methods,
                                                   conda.bin = NULL,
                                                   conda.env.name = NULL,
                                                   conda.env.path = NULL) {
  method <- tolower(method)
  method <- match.arg(method)
  slot(object = x, name = method) <- updateEnv(slot(object = x, name = method),
                                               conda.bin = conda.bin,
                                               conda.env.name = conda.env.name,
                                               conda.env.path = conda.env.path)
  x
})

setMethod("show", "CondaEnvSlot", function(object) {
  cat(c("An invalid", "A valid")[isValid(object) + 1],
      " element of a conda environement (", sQuote(object@value), ")\n",
      sep = "")
})
setMethod("show", "CondaEnv", function(object) {
  # header1 <- cli::bg_blue(cli::ansi_align(cli::col_br_white(cli::style_bold("Conda environment")), align="center"))
  header2 <- cli::rule(left = cli::style_bold("Conda environment"),
                       right = sprintf("%s (v%s)", "SeuratIntegrate", getNamespaceVersion(.packageName)[[1]]),
                       line_col = "gray60", col= "royalblue")
  good_tick <- cli::col_green(cli::symbol$tick)
  bad_cross <- cli::col_red(cli::symbol$cross)
  sep <- paste(" ", cli::style_strikethrough("    "), "  ")
  is.valid <- good_tick
  env.name <- cli::col_blue(cli::style_dim("R-based method"))
  current.status <- paste0(is.valid, "  ", object$method, sep, env.name,
                           "\n", collapse = " ")
  if (object$needs.conda) {
    is.valid <- ifelse(isValid(object), good_tick, bad_cross)
    env.name <- ifelse(isValid(object), cli::col_blue(object$conda.env.name$value), bad_cross)
    env.path <- ifelse(isValid(object), cli::col_blue(dirname(object$conda.env.path$value)), bad_cross)
    current.status <- paste0(is.valid, "  ", object$method, sep, env.name,
                             sep, env.path, "\n", collapse = " ")
  }
  cat(header2, "\n ", current.status, sep = "")
})
setMethod("show", "CondaEnvManager", function(object) {
  good_tick <- cli::col_green(cli::symbol$tick)
  bad_cross <- cli::col_red(cli::symbol$cross)
  is.valid <- ifelse(isValid(object), good_tick, bad_cross)
  methods <- cli::style_bold(names(is.valid))
  needs.conda <- sapply(object, function(x) x@needs.conda)[names(is.valid)]

  env.name <- sapply(object, function(x) x@conda.env.name@value)[names(is.valid)]
  env.name[!needs.conda] <-  cli::style_dim("R-based method")
  env.name <- cli::col_blue(env.name)
  env.name[!isValid(object)] <- bad_cross

  env.path <- sapply(object, function(x) x@conda.env.path@value)[names(is.valid)]
  env.path[!needs.conda] <-  cli::style_dim("R-based method")
  env.path <- cli::col_grey(env.path)
  env.path[!isValid(object)] <- bad_cross
  env.path[isValid(object) & needs.conda] <- dirname(env.path[isValid(object) & needs.conda])
  sep <- paste(" ", cli::style_strikethrough("    "), "  ")
  current.status <- paste0(is.valid, "  ",
                           cli::ansi_align(methods, max(cli::ansi_nchar(methods))), sep,
                           cli::ansi_align(env.name, max(cli::ansi_nchar(env.name))),
                           sep, env.path, "\n", collapse = " ")
  # header1 <- cli::rule(center =  cli::style_bold("Conda environments manager"),
  #                      line = " ", line_col = "black", col = "gray30")
  header1 <- cli::bg_blue(cli::ansi_align(cli::col_br_white(cli::style_bold("Conda environments manager")), align="center"))
  header2 <- cli::rule(left = cli::style_bold("Integration methods"),
                       right = sprintf("%s (v%s)", "SeuratIntegrate", getNamespaceVersion(.packageName)[[1]]),
                       line_col = "gray60", col= "royalblue")
  cat(header1, "\n", header2, "\n\n ", current.status, sep = "")
})


################################################################################
##########################                            ##########################
#######                        set R base methods                        #######
###                                                                          ###
#                                                                              #

#' @importFrom methods slotNames
#' @export
as.list.CondaEnvSlot <- function(x, ...) {
  sapply(slotNames(x), slot, object = x, simplify = FALSE, USE.NAMES = TRUE)
}
#' @export
as.list.CondaEnv <- function(x, ...) {
  unlist(sapply(
    sapply(slotNames(x), slot, object = x, simplify = FALSE, USE.NAMES = TRUE),
    as.list), recursive = FALSE)
}
#' @export
as.list.CondaEnvManager <- function(x, recursive = FALSE, ...) {
  l <- sapply(slotNames(x), slot, object = x, simplify = FALSE, USE.NAMES = TRUE)
  if (recursive) {
    l <- lapply(l, as.list)
  }
  l
}
#' @export
as.data.frame.CondaEnvSlot <- function(x, ...) {
  as.data.frame(as.list(x), ...)
}
#' @export
as.data.frame.CondaEnv <- function(x, ...) {
  as.data.frame(as.list(x), row.names = x@method)
}
#' @export
as.data.frame.CondaEnvManager <- function(x, ...) {
  do.call(rbind, lapply(x, as.data.frame))
}


#' @importFrom utils .DollarNames
#' @export
.DollarNames.CondaEnvSlot <- function(x, pattern="") {
  return("value")
}
#' @importFrom utils .DollarNames
#' @export
.DollarNames.CondaEnv <- function(x, pattern="") {
  return(c("method", "needs.conda", "conda.bin", "conda.env.name",
           "conda.env.path", "is.valid"))
}
#' @importFrom utils .DollarNames
#' @export
.DollarNames.CondaEnvManager <- function(x, pattern="") {
  return(known.methods)
}

#' @export
"$.CondaEnvSlot" <- function(object, name = NULL) {
  name <- name %||% "value"
  res <- NULL
  if (name == "value") {
    res <- slot(object = object, name = "value")
  }
  return(res)
}
#' @export
"[[.CondaEnvSlot" <- `$.CondaEnvSlot`

#' @export
"$.CondaEnv" <- function(object, name = NULL) {
  name <- intersect(name, c("method", "needs.conda", "conda.bin",
                            "conda.env.name", "conda.env.path"))
  if (length(name) == 1) {
    res <- slot(object = object, name = name)
  } else if (length(name) == 0) {
    warning("Unknown", call. = F, immediate. = T)
    res <- NULL
  } else {
    res <- object[[name]]
  }
  return(res)
}
#' @export
"[[.CondaEnv" <- `$.CondaEnv`

#' @export
"$.CondaEnvManager" <- function(object, name) {
  name <- intersect(name, known.methods)
  if (length(name) == 1) {
    res <- slot(object = object, name = name)
  } else if (length(name) == 0) {
    warning("Unknown method", call. = F, immediate. = T)
    res <- NULL
  } else {
    res <- object[name]
  }
  return(res)
}
#' @export
"[[.CondaEnvManager" <- `$.CondaEnvManager`

#' @export
"[.CondaEnvManager" <- function(object, name = NULL) {
  name <- name %||% known.methods
  return(as.data.frame(object)[intersect(name, known.methods),])
}

#' @export
"$<-.CondaEnvSlot" <- function(object, name, value) {
  if (name == "value") {
    slot(object = object, name = "value") <- value
    slot(object = object, name = "valid") <- FALSE
  } else {
    warning("Only ", sQuote("value"), " can be updated. Object unchanged",
            call. = F, immediate. = T)
  }
  return(object)
}
#' @importFrom methods is
#' @export
"$<-.CondaEnvManager" <- function(object, name, value) {
  name <- intersect(name, known.methods)
  if (length(name) == 0) {
    warning("Unknown method. Object unchanged", call. = F, immediate. = T)
  } else if (length(name) > 1) {
    warning("Can only modify 1 method. Object unchanged", call. = F, immediate. = T)
  } else {
    value <- value %||% .blankCondaManager()[[name]]
    if (! is(value, "CondaEnv")) {
      warning(sQuote(class(value)), " is not ", sQuote("CondaEnv"),
              ". Object unchanged", call. = F, immediate. = T)
    } else if (name != value@method) {
      warning("Discrepancy between ", sQuote(name), " and provided value's method (",
              sQuote(value@method),"). Object unchanged", call. = F, immediate. = T)
    } else {
      slot(object = object, name = name) <- value
    }
  }
  return(object)
}
#' @export
"[<-.CondaEnvManager" <- `$<-.CondaEnvManager`
#' @export
"[[<-.CondaEnvManager" <- `$<-.CondaEnvManager`
