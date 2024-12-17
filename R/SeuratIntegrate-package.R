#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang %||% check_installed abort
#' @importFrom SeuratObject %iff%
#' @import rlang
## usethis namespace: end

#' @name integration-method
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' (or an \code{\link[SeuratObject]{Assay5}} object if
#' not called by \code{\link[Seurat]{IntegrateLayers}})
#' @param orig \code{\link[SeuratObject]{DimReduc}} object. Not to be set
#' directly when called with \code{\link[Seurat]{IntegrateLayers}}, use
#' \code{orig.reduction} argument instead
#' @param groups A \bold{named} data frame with grouping information. Preferably
#' one-column when \code{groups.name = NULL}
#' @param groups.name Column name from \code{groups} data frame that stores
#' grouping information. If \code{groups.name = NULL}, the first column is used
#' @param layers Name of the layers to use in the integration
#' @param scale.layer Name of the scaled layer in \code{Assay}
#' @param seed.use An integer to generate reproducible outputs.
#' Set \code{seed.use = NULL} to disable
#' @param verbose Print messages. Set to \code{FALSE} to disable
#' @param features Vector of feature names to input to the integration method.
#' When \code{features = NULL} (default), the
#' \code{\link[SeuratObject]{VariableFeatures}} are used. To pass all features,
#' use the output of \code{\link[SeuratObject]{Features}()}
#' @param new.reduction Name of the new integrated dimensional reduction
#' @param reduction.key Key for the new integrated dimensional reduction
NULL
