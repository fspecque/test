#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang %||%
#' @importFrom SeuratObject %iff%
## usethis namespace: end
#' @name integration-method
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' (or an \code{\link[SeuratObject]{Assay5}} object if
#' not called by \code{\link[Seurat]{IntegrateLayers}})
#' @param orig \code{\link[SeuratObject]{DimReduc}} object. Not to be set
#' directly when called with \code{\link[Seurat]{IntegrateLayers}}, use
#' \code{orig.reduction} argument instead
#' @param groups A data frame with grouping information. Should be one-column
#' when \code{groups.name = NULL}
#' @param groups.name Column name from \code{groups} data frame that stores
#' grouping information
#' @param layers Name of layers to use in the integration
#' @param verbose Print messages. Set to \code{FALSE} to disable
NULL
