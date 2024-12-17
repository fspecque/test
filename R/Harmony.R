# Fixes:
## choice of number of dimensions used
## choice of number of threads    used
#' Run Harmony on Seurat's \link[SeuratObject]{Assay5} object through \code{\link[Seurat]{IntegrateLayers}}
#' @description
#' A wrapper to run \code{\link[harmony:RunHarmony]{Harmony}} on multi-layered Seurat V5 object
#'
#' Can be called via \code{SeuratIntegrate::HarmonyIntegration()} or
#' \code{HarmonyIntegration.fix()}
#'
#' @inheritParams integration-method
#' @inheritParams harmony::RunHarmony.default
#' @param dims Dimensions of dimensional reduction to use for integration.
#' All used by default
#' @param layers Ignored unless \code{groups = NULL}, then used to create
#' grouping variable to correct batch-effect.
#' @param scale.layer Ignored
#' @param features Ignored
#' @param key Prefix for the dimension names computed by harmony.
#' @param ... Ignored for \code{HarmonyIntegration()}, or all of the above for
#' \code{HarmonyIntegration.fix()}
#'
#' @return The function itself returns a list containing:
#' \itemize{
#'   \item a new DimReduc of name \code{reduction.name} (key set to
#'   \code{reduction.key}) with corrected cell embeddings matrix of
#'   \code{length(dims)} columns.
#' }
#' When called via \code{\link[Seurat]{IntegrateLayers}}, a Seurat object with
#' the new reduction is returned
#'
#' @importFrom SeuratObject Embeddings CreateDimReducObject DefaultAssay
#' @importFrom harmony harmony_options RunHarmony
#' @importFrom stats sd
#'
#' @export
#' @note This function requires the
#' \href{https://portals.broadinstitute.org/harmony/}{\pkg{harmony}} package
#' to be installed
#'
#' @references Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F.,
#' Wei, K., Baglaenko, Y., Brenner, M., Loh, P. & Raychaudhuri, S. Fast,
#' sensitive and accurate integration of single-cell data with Harmony.
#' Nat Methods 16, 1289–1296 (2019). \href{https://doi.org/10.1038/s41592-019-0619-0}{DOI}
#'
#' @examples
#' \dontrun{
#' # Preprocessing
#' obj <- UpdateSeuratObject(SeuratData::LoadData("pbmcsca"))
#' obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj)
#'
#' # After preprocessing, we integrate layers based on the "Method" variable:
#' obj <- IntegrateLayers(object = obj, method = SeuratIntegrate::HarmonyIntegration,
#'                        verbose = TRUE)
#'
#' # We can also change parameters such as the batch-effect variable.
#' # Here we change the groups variable, the number of dimension used from the original
#' # PCA and minor options from `harmony_options()`:
#' harmonyOptions <- harmony::harmony_options()
#' harmonyOptions$max.iter.cluster <- 10   #  20 by default
#' harmonyOptions$block.size <- .1         # .05 by default
#' obj <- IntegrateLayers(object = obj, method = SeuratIntegrate::HarmonyIntegration,
#'                        dims = 1:30, plot_convergence = TRUE,
#'                        groups = obj[[]]$Experiment,
#'                        new.reduction = "harmony_custom",
#'                        .options = harmonyOptions, verbose = TRUE)
#' }
HarmonyIntegration <- function(
    object,
    orig,
    groups = NULL,
    groups.name = NULL,
    layers = NULL,
    scale.layer = 'scale.data',
    features = NULL,
    new.reduction = 'harmony',
    dims = NULL,     # harmonize with CCA/rPCA
    key = 'harmony_',
    seed.use = 42L,
    theta = NULL,
    sigma = 0.1,
    lambda = NULL,
    nclust = NULL,
    ncores = 1L,
    max_iter = 10,
    early_stop = TRUE,
    plot_convergence = FALSE,
    # return_object = FALSE,
    .options = harmony_options(),
    verbose = TRUE,
    ...
) {
  check_installed(
    pkg = "harmony",
    reason = "for running integration with Harmony"
  )
  if (!inherits(x = object, what = c('StdAssay', 'SCTAssay'))) {
    abort(message = "'object' must be a v5 or SCT assay")
  } else if (!inherits(x = orig, what = 'DimReduc')) {
    abort(message = "'orig' must be a dimensional reduction")
  }

  data_mat <- Embeddings(object = orig)
  dims <- dims %||% 1:ncol(data_mat)
  if (! all(dims %in% 1:ncol(data_mat))) {
    warning(sprintf(
      "'dims' must not be greater than number of dimensions in 'orig'. Setting to 1:%d",
      ncol(data_mat)), call. = T)
    dims <- 1:ncol(data_mat)
  }

  layers <- Layers(object = object, search = layers %||% 'counts')
  groups <- groups %||% CreateIntegrationGroups(object = object,
                                                layers = layers,
                                                scale.layer = scale.layer)
  if (! inherits(x = groups, what = "data.frame")) {
    # groups is supposedly a vector, a matrix or a list
    groups <- as.data.frame(groups)
  }
  groups.name <- groups.name %||% colnames(groups)[1]
  groups.name <- intersect(colnames(groups), groups.name)
  if (! length(x = groups.name)) {
    abort(message = "'groups.name' not in 'groups' data frame")
  }
  if (length(x = groups.name) > 1) {
    warning(paste("more 'groups.name' that expected. Using the first one",
                  sQuote(x = groups.name[1])), call. = FALSE, immediate. = TRUE)
    groups.name <- groups.name[1]
  }
  msg <- sprintf("Using groups in %s encompassing %d levels (%s)\n",
                 sQuote(x = groups.name),
                 length(unique(groups[,groups.name, drop=T])),
                 paste0(sQuote(unique(groups[,groups.name, drop=T])),
                        collapse = ", "))
  message(msg[verbose], appendLF = FALSE)

  # Run Harmony
  seed.use <- ifelse(isFALSE(seed.use), NULL, seed.use)
  seed.use %iff% set.seed(seed.use)
  harmony.embed <- RunHarmony(
    data_mat = data_mat[,dims,drop=F],
    meta_data = groups,
    vars_use = groups.name,
    theta = theta,
    sigma = sigma,
    lambda = lambda,
    nclust = nclust,
    max_iter = max_iter,
    early_stop = early_stop,
    ncores = ncores,
    plot_convergence = plot_convergence,
    return_object = FALSE,
    verbose = verbose,
    .options = .options
  )
  rownames(x = harmony.embed) <- Cells(x = orig)
  # TODO add feature loadings from PCA
  dr <- suppressWarnings(expr = CreateDimReducObject(
    embeddings = harmony.embed,
    stdev = as.numeric(apply(harmony.embed, 2, sd)),
    key = key,
    assay = DefaultAssay(object = orig)
  ))
  output.list <- list(dr)
  names(output.list) <- c(new.reduction)
  return(output.list)
}

attr(x = HarmonyIntegration, which = 'Seurat.method') <- 'integration'

#' @rdname HarmonyIntegration
#' @export
HarmonyIntegration.fix <- function(...) {
  HarmonyIntegration(...)
}
attr(x = HarmonyIntegration.fix, which = 'Seurat.method') <- 'integration'
