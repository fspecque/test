# define 'numeric_lisi' -> numeric vector with @N slot storing the number of
# categories in the variable used to compute lisi score
# => used for scaling after saving to Misc

#' `numeric_lisi` is built to expand `numeric` with a sticky integer
#' @description
#' A numeric vector with `@N` slot storing the number of categories in the
#' variable used to compute the lisi score. Used for scaling after saving to Misc
#'
#' @keywords internal
#' @noRd
setClass('numeric_lisi', contains = "numeric",
         slots = c(N = "integer"),
         prototype = structure(numeric(), N = NA_integer_))


#' `numeric_lisi` constructor
#' @keywords internal
#' @noRd
numeric_lisi <- function(x, N = NA_integer_) {
  if (is.numeric(N)) {
    N <- as.integer(N)
  } else if (suppressWarnings(!is.na((N_ <- as.integer(N))))) {
    N <- N_
  } else { N <- NA_integer_}
  new('numeric_lisi', x, N = N)
}


#' Get and set `@N` slot of a `numeric_lisi`
#' @keywords internal
#' @noRd
setGeneric("N", function(x) standardGeneric("N"))
setGeneric("N<-", function(x, value) standardGeneric("N<-"))

setMethod("N", "numeric_lisi", function(x) x@N)
setMethod("N<-", "numeric_lisi", function(x, value) {
  if (is.numeric(value)) {
    value <- as.integer(value)
  } else if (suppressWarnings(!is.na((value_ <- as.integer(value))))) {
    value <- value_
  }
  x@N <- value
  validObject(x)
  x
})

#' Check validity of `numeric_lisi` (`@N` slot)
#' @keywords internal
#' @noRd
setValidity('numeric_lisi', function(object) {
  if (length(object@N) > 1) {
    "number of categories (@N) should be a single integer"
  } else {
    TRUE
  }
})


#' Overload print method for `numeric_lisi` (hide `@N` slot)
#' @keywords internal
#' @noRd
setMethod('show', 'numeric_lisi', function(object) {
  attributes(object) <- NULL
  print(unclass(object))
})


# Coerce types
#' Coercing from and to `numeric_lisi` (hide `@N` slot)
#' @keywords internal
#' @noRd

##  use `as('from', 'to')`
setAs('numeric', 'numeric_lisi', function(from) numeric_lisi(from))
setAs('double', 'numeric_lisi', function(from) numeric_lisi(from))
setAs('ANY', 'numeric_lisi', function(from) numeric_lisi(as.numeric(from)))
setAs('numeric_lisi', 'numeric', function(from) as.numeric(unclass(from)))

##  use `as.to(from)`
setGeneric("as.numeric_lisi", function(x) standardGeneric("as.numeric_lisi"))
setMethod('as.numeric_lisi', 'numeric', function(x) {
  numeric_lisi(x)
})
setMethod('as.numeric_lisi', 'numeric_lisi', identity)
setMethod('as.numeric_lisi', 'list', function(x) {
  as.numeric_lisi(as.numeric(x))
})
setMethod('as.numeric_lisi', 'logical', function(x) {
  as.numeric_lisi(as.numeric(x))
})
setMethod('as.numeric', 'numeric_lisi', function(x) {
  as.numeric(unclass(x))
})


#' Indexing method (keep `@N` slot)
#' @keywords internal
#' @noRd
setMethod('[', 'numeric_lisi', function(x, i) {
  numeric_lisi(callNextMethod(), N(x))
})


#' Replicating method (keep `@N` slot)
#' @keywords internal
#' @noRd
setMethod('rep', 'numeric_lisi', function(x, ...) {
  numeric_lisi(callNextMethod(), N(x))
})


#' Helper method to pick `@N` slot and check compatibility between `numeric_lisi`
#' @keywords internal
#' @noRd
pick_N <- function(x, y) {
  picked.N <- N(x)
  if (is.na(N(x))) {
    picked.N <- N(y)
  } else if (!is.na(N(y))) {
    if (N(x) != N(y)) {
      abort('incompatible number of categories between lisi result vectors')
    }
  }
  return(picked.N)
}

#' Value replacement method (check `@N` slot compatibility)
#' @keywords internal
#' @noRd
setMethod('[<-', 'numeric_lisi', function(x, i, value) {
  if ('numeric_lisi' %in% class(value)) {
    picked_N <- pick_N(x, value)
  } else { picked_N <- N(x)}
  numeric_lisi(callNextMethod(), picked_N)
})


#' Casting and coercion support
#'
#'  @description
#'  Double dispatch methods to support [vctrs::vec_cast()] and [vctrs::vec_ptype2].
#'  Used in [tibble::add_row()] inside `AddMiscIntegrations()` chiefly
#'
#' @inheritParams vctrs::vec_cast
#'
#' @method vec_cast numeric_lisi
#' @exportS3Method vctrs::vec_cast
#' @importFrom vctrs vec_cast vec_ptype2
#' @keywords internal
#' @noRd
vec_cast.numeric_lisi <- function(x, to, ...) UseMethod("vec_cast.numeric_lisi")

#' @keywords internal
#' @noRd
#' @method vec_ptype2 numeric_lisi
#' @exportS3Method vctrs::vec_ptype2
vec_ptype2.numeric_lisi <- function(x, y, ...) UseMethod("vec_ptype2.numeric_lisi")

#' @method vec_ptype2.numeric_lisi numeric
#' @export
vec_ptype2.numeric_lisi.numeric <- function(x, y, ...) {
  numeric_lisi(numeric(), N = N(x))
}
#' @method vec_ptype2.numeric_lisi double
#' @export
vec_ptype2.numeric_lisi.double <- function(x, y, ...) {
  numeric_lisi(numeric(), N = N(x))
}

#' @method vec_ptype2.numeric_lisi numeric_lisi
#' @export
vec_ptype2.numeric_lisi.numeric_lisi <- function(x, y, ...) {
  numeric_lisi(numeric(), N = pick_N(x, y))
}

#' @method vec_cast.numeric_lisi numeric_lisi
#' @export
vec_cast.numeric_lisi.numeric_lisi <- function(x, to, ...) {
  numeric_lisi(x, N = pick_N(x, to))
}
