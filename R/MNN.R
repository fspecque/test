#' Run classical MNN on Seurat's \link[SeuratObject]{Assay5} object through \code{\link[Seurat]{IntegrateLayers}}
#'
#' A wrapper to run \code{\link[batchelor]{mnnCorrect}} on multi-layered Seurat V5 object
#'
#' @inheritParams integration-method
#' @param groups Ignored
#' @param features Either a list of features to use when calculating batch
#' correction, or a number (2000 by default) of variable features to select.
# @param reduction.name Name to store resulting DimReduc object
# @param reduction.key Key for resulting DimReduc
#' @param reconstructed.assay Name for the assay containing the low-rank
#' reconstruction of the expression matrix.
#' @param ... Extra parameters passed to \code{\link[batchelor]{mnnCorrect}}
#'
#' @return A Seurat object merged from the objects in \code{object.list} and a
#' new DimReduc of name \code{reduction.name} (key set to \code{reduction.key})
#' with corrected embeddings matrix as well as the rotation matrix used for the
#' PCA stored in the feature loadings slot. Also returns an expression matrix
#' reconstructed from the low-rank approximation in the
#' \code{reconstructed.assay} assay; all other metadata info
#' \code{\link[batchelor]{mnnCorrect}} is stored in the \code{tool} slot,
#' accessible with \code{\link[Seurat]{Tool}}
#'
#' @importFrom SeuratObject Layers LayerData CreateAssayObject
#' @importFrom Seurat DefaultAssay DefaultAssay<- SelectIntegrationFeatures VariableFeatures VariableFeatures<-
#' as.SingleCellExperiment CreateDimReducObject Tool<- LogSeuratCommand CreateSeuratObject SelectIntegrationFeatures5
#' @importFrom batchelor mnnCorrect
#' @importFrom SummarizedExperiment assay
#'
#' @export
#' @note This function requires the
#' \href{https://bioconductor.org/packages/release/bioc/html/batchelor.html}{\pkg{batchelor}} package
#' to be installed
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
#' # After preprocessing, we integrate layers:
#' obj <- IntegrateLayers(object = obj, method = MNNIntegration,
#'   new.reduction = 'integrated.mnn', verbose = FALSE)
#'
#' # We can also add parameters specific to mnnCorrect.
#' # Here we set `k` to specify the number of nearest neighbors
#' # to use when identifying MNNs:
#' obj <- IntegrateLayers(object = obj, method = MNNIntegration,
#'   new.reduction = 'integrated.mnn', k = 15, verbose = FALSE)
#' }
#'
#' @seealso \code{\link[batchelor]{mnnCorrect}} \code{\link[Seurat]{Tool}}
#'
### Almost copy-past from FastMNNIntegration
MNNIntegration <- function(
    object,
    orig = NULL,
    groups = NULL,
    layers = NULL,
    scale.layer = NULL,
    features = 2000,
    # reduction.key = "mnn_",
    reconstructed.assay = "mnn.reconstructed",
    verbose = TRUE,
    ...
) {
  check_installed(
    pkg = "batchelor",
    reason = "for running integration with mnnCorrect"
  )

  varargs <- list(...)

  object <- CreateSeuratObject(object)
  if (is.numeric(x = features)) {

    message(paste("Computing", features, "integration features\n")[verbose],
            appendLF = F)

    features <- SelectIntegrationFeatures5(object = object, features = features)
  }
  layers <- Layers(object, search = layers %||% 'data')
  message("Converting layers to SingleCellExperiment\n"[verbose], appendLF = F)

  objects.sce <- lapply(
    X = layers,
    FUN = function(x, f) {
      return(as.SingleCellExperiment(
        x = subset(x = object,
                   features = f,
                   cells = colnames(LayerData(object, layer = x)))
      )
      )
    },
    f = features
  )

  message("Running MNN correction\n"[verbose], appendLF = F)
  out <- do.call(
    what = mnnCorrect,
    args = c(
      objects.sce,
      varargs[intersect(names(varargs), args$combat$mnnCorrect)]
    )
  )
  message("Done MNN correction\n"[verbose], appendLF = F)

  # Add reconstructed matrix (gene x cell)
  reconstructed_assay <- CreateAssayObject(
    data = choose_matrix_format(mat = assay(x = out)),
  )
  # Add variable features
  VariableFeatures(object = reconstructed_assay) <- features
  #Tool(object = object) <- S4Vectors::metadata(x = out)
  #object <- LogSeuratCommand(object = object)
  output.list <- setNames(list(reconstructed_assay), reconstructed.assay)
  return(output.list)
}

attr(x = MNNIntegration, which = 'Seurat.method') <- 'integration'
