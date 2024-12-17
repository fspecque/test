#' Integrate layers using one or multiple integration method(s)
#' @description
#' Integrate layers of a Seurat object using one or more integration methods.
#'
#' Available integration methods are listed at the bottom of this page.
#' \code{DoIntegrate()} works best with \pkg{SeuratIntegrate}'s methods.
#'
#' @inheritSection Seurat::IntegrateLayers Integration Method Functions
#'
#' @param object A Seurat object
#' @param ... one or more integration method \strong{call}(s) to correct batch
#' effects with. Must be of the form \code{package::MethodIntegration()}. It is
#' recommended to use \code{::} because it is safer in case of namespace
#' collision or if the package is not attached nor loaded. \strong{Don't forget
#' the parentheses}.
#' @param use.hvg whether to use highly variable genes. \code{FALSE} causes all
#' the features present in the assay to be used.
#' @param use.future whether to use \pkg{future} to run integrations in
#' a background session. Useful when python-based algorithms are invoked.
#' @param future.globals.size maximum allowed size (in bytes) of global
#' variables to export. By default, uses the value of the option
#' "future.globals.maxSize". If it is \code{NULL}, thrice the size of the Seurat
#' object is used. Inoperative when \code{use.future = FALSE}
#'
#' @export
#'
#' @return the updated Seurat object enriched with the integration methods'
#' outputs.
#'
#' @importFrom rlang enquos as_label
#'
#' @details
#' Each call to an integration method \strong{require parentheses}. Parameter
#' values specific to each method can be enclosed between these brackets,
#' although the defaults arguments should be good enough for most cases. Note
#' that for each integration method, \strong{the argument values specified in
#' its call supersede \code{DoIntegrate}'s internal computations}. For instance,
#'
#' \preformatted{
#' DoIntegrate(seu,
#'             SeuratIntegrate::CombatIntegration(features = Features(seu)),
#'             Seurat::CCAIntegration(),
#'             use.hvg = TRUE)
#' }
#'
#' forces ComBat but not CCA to use all features instead of the variable ones.
#'
#'
#' @note The desired value of parameters \code{use.hvg}, \code{use.future} and
#' \code{future.globals.size} can be different for each method called. Hence,
#' they accept vectors with more than one element. They must be in the same
#' order as the integration method calls.
#' @note This unconventional way of calling methods presents a few advantages:
#' With a single call to \code{DoIntegrate}, you can perform multiple
#' integrations at once, while preserving the flexibility of an individual
#' method call by fine-tuning integration-specific parameters.
#'
DoIntegrate <- function (object, ..., use.hvg = TRUE, use.future = TRUE,
                         future.globals.size = getOption("future.globals.maxSize")) {
  integrations <- rlang::enquos(..., .ignore_empty = 'all', .ignore_null = 'all')
  n.integrations <- length(integrations)
  use.hvg <- rep(use.hvg, length.out = n.integrations)
  use.future <- rep(use.future, length.out = n.integrations)
  future.globals.size <- future.globals.size %iff% rep(future.globals.size, length.out = n.integrations)

  for (int in 1:n.integrations) {
    cat(sprintf('Integration %d in %d: ', int, n.integrations))
    tryCatch({
    object <- DoIntegrateSingle(
      object = object, method = integrations[[int]], use.hvg = use.hvg[int],
      use.future = use.future[int],
      future.globals.size = future.globals.size %iff% future.globals.size[int])
    }, error = function(e) {
      warning("`", sub("\\(.*", "", as_label(integrations[[int]])), "`",
              ": failed, because ", dQuote(e$message),
              "\n -> Skipping this integration",
              call. = FALSE, immediate. = TRUE)
    })
  }

  return(object)
}

#' @importFrom rlang quo_get_expr call_ns call_name call2 call_match call_args eval_tidy maybe_missing fn_fmls_names call_modify quo quo_set_expr
#' @importFrom SeuratObject DefaultAssay Layers Features VariableFeatures DefaultDimReduc
#' @importFrom Seurat SelectSCTIntegrationFeatures
#' @importFrom parallelly supportsMulticore
#' @importFrom future sequential multisession multicore plan %<-% %seed% %packages%
#' @importFrom utils object.size
#' @keywords internal
#' @noRd
DoIntegrateSingle <- function (object, method, use.hvg = TRUE, use.future = TRUE,
                         future.globals.size = getOption("future.globals.maxSize"))
{
  method_expr <- quo_get_expr(method)
  is_a_call <- tryCatch(call_name(method_expr) %iff% TRUE %||% FALSE,
                        error = function(e) FALSE)
  method_expr <- if (! is_a_call)  call2(method_expr) else method_expr
  method_fun <- paste0(c(call_ns(method_expr),
                         call_name(method_expr)),
                       collapse = '::')
  cat(sprintf('integrating using %s\n', sQuote(method_fun)))
  method_env <- eval(parse(text = method_fun))
  method_expr <- call_match(method_expr, method_env, defaults = FALSE)

  method_args <- call_args(method_expr)#lapply(call_args(method_expr), quo)
  if (use.future && ! supportsMulticore()) {
    # suboptimal but fixes object no found for future plan(multisession)
    method_args <- lapply(method_args, eval_tidy)
  }
  assay <- method_args$assay %||% DefaultAssay(object = object)
  layers <- eval_tidy(method_args$layers)
  scale.layer <- eval_tidy(method_args$scale.layer) %||% 'scale.data'
  scale.layer <- Layers(object = object, search = scale.layer, assay = assay)
  orig.use <- maybe_missing(eval_tidy(method_args$orig), default = NULL)
  features <- eval_tidy(method_args$features)

  features <- features %||% if (!use.hvg) {
    Features(object[[assay]], layer = layers)
  } else if (inherits(x = object[[assay]], what = "SCTAssay")) {
    mean.features <- mean(
      sapply(levels(object[['SCT']]), function(lvl) {
        length(VariableFeatures(object[['SCT']], layer = lvl))
      }, simplify = TRUE, USE.NAMES = FALSE)
    )
    SelectSCTIntegrationFeatures(object = object, assay = assay,
                                 nfeatures = mean.features)
  } else if (inherits(x = object[[assay]], what = "StdAssay")) {
    VariableFeatures(object = object, assay = assay)
  } else {
    abort(message = "'assay' must be a v5 or SCT assay")
  }

  features <- intersect(features, Features(object[[assay]], layer = layers))
  if (!length(x = features)) {
    abort(message = "None of the features provided are found in this assay")
  }
  if ('orig' %in% fn_fmls_names(method_env)) {
    orig.name <- NULL
    if (orig.use %iff% is.character(orig.use) %||% TRUE) {
      orig.name <- orig.use %||% DefaultDimReduc(object = object, assay = assay)
      orig.use <- object[[orig.name]]
    }
  }
  assay.use <- object[[assay]]
  method_env_name <- tryCatch({
    call_ns(method_expr) %||% environmentName(environment(fun = get(method_fun)))},
    error = function(e) "Unknown")
  if (inherits(assay.use, 'SCTAssay') && method_env_name == .packageName) {
    groups <- CreateIntegrationGroups(assay.use)
    assay.use <- suppressMessages(suppressWarnings(
      split(assay.use, f = groups$group, layers = c("counts", "data"))))
  }
  method_expr <- call_modify(method_expr, !!! method_args)
  method_expr <- call_modify(method_expr, object = quo(assay.use), assay = assay,
                             orig = quo(orig.use), layers = layers,
                             scale.layer = scale.layer, features = quo(features))
  if (use.future) {
    future.globals.size <- future.globals.size %||% unclass(object.size(object)) * 3
    oopts <- options(future.globals.maxSize = future.globals.size)  ## 1.0 GB
    oplan <- plan()
    on.exit(options(oopts))
    on.exit(plan(oplan))

    future.strat <- if (supportsMulticore()) multicore else multisession
    plan(future.strat)
    value %<-% {
      eval_tidy(quo_set_expr(method, method_expr))
    } %seed% 42L %packages% unique(c('Seurat', 'SeuratIntegrate', call_ns(method_expr)))
  } else {
    value <- eval_tidy(quo_set_expr(method, method_expr))
  }

  for (i in names(x = value)) {
    object[[i]] <- value[[i]]
  }
  return(object)
}
