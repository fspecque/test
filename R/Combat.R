#' Run ComBat on Seurat's \link[SeuratObject]{Assay5} object through \code{\link[Seurat]{IntegrateLayers}}
#' @description
#' A wrapper to run \code{\link[sva:ComBat]{ComBat}} or
#' \code{\link[sva:ComBat_seq]{ComBat_seq}} on multi-layered Seurat V5 object
#'
#'
#' @inheritParams integration-method
#' @param reconstructed.assay Name for the \code{assay} containing the corrected
#' expression matrix
#' @param key.assay Optional key for the new combat assay. Format: "[:alnum:]*_"
#' @param combat.function ComBat implementation to use. One of
#' \link[sva:ComBat]{combat}, \link[sva:ComBat_seq]{combat_seq}. Note that
#' ComBat_seq is an improved model from ComBat but requires a dense matrix.
#' Sparse to dense matrix conversion can be memory-intensive.
#' @param use.scaled By default the layer passed to the \code{layer} argument is
#' used. When \code{use.scaled = TRUE}, the \code{scale.layer} is input to ComBat.
#' @param ... Additional arguments passed on to \link[sva:ComBat]{ComBat} or
#' \link[sva:ComBat_seq]{ComBat_seq}.
#'
#' @return The function itself returns a list containing:
#' \itemize{
#'   \item a new Assay of name \code{reconstructed.assay} (key set to
#'   \code{assay.key}) with corrected cell counts.
#' }
#' When called via \code{\link[Seurat]{IntegrateLayers}}, a Seurat object with
#' the new assay is returned
#'
#' @importFrom SeuratObject Layers JoinLayers GetAssayData CreateAssayObject
#' SetAssayData
#' @importFrom sva ComBat ComBat_seq
#'
#' @export
#' @note This function requires the
#' \href{https://bioconductor.org/packages/release/bioc/html/sva.html}{\pkg{sva}
#' (Surrogate Variable Analysis)} package to be installed
#'
#' @references Johnson, W. E., Li, C. & Rabinovic, A. Adjusting batch effects in
#' microarray expression data using empirical Bayes methods. Biostatistics 8,
#' 118–127 (2006). \href{https://doi.org/10.1093/biostatistics/kxj037}{DOI}
#'
#' Zhang, Y., Parmigiani, G. & Johnson, W. E. ComBat-seq: batch effect
#' adjustment for RNA-seq count data. NAR Genomics and Bioinformatics 2 (2020).
#' \href{https://doi.org/10.1093/nargab/lqaa078}{DOI}
#'
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
#' obj <- IntegrateLayers(object = obj, method = CombatIntegration,
#'                        verbose = TRUE, layers = "data", scale.layer = NULL,
#'                        features = VariableFeatures(
#'                          FindVariableFeatures(obj, nfeatures = 5e3)
#'                        ))
#'
#' # We can also change parameters such as the input data.
#' # Here we use the scale data, the ComBat implementation and we use the cell
#' # labels as a "biological condition of interest" (/!\ long):
#'
#' obj <- IntegrateLayers(object = obj,  method = CombatIntegration,
#'                        verbose = TRUE, features = VariableFeatures(obj),
#'                        use.scaled = FALSE, combat.function = 'combat_seq',
#'                        group = obj[[]]$CellType, groups = obj[[]],
#'                        groups.name = "Method", layers = "counts")
#' }

CombatIntegration <- function(
    object,
    orig = NULL,
    groups = NULL,
    groups.name = NULL,
    layers = 'data',
    scale.layer = 'scale.data',
    features = NULL,
    reconstructed.assay = "combat.reconstructed",
    key.assay = "combat_",
    combat.function = c("combat", "combat_seq"),
    use.scaled = FALSE,
    verbose = TRUE,
    ...
) {
  combat.function <- tolower(combat.function)
  combat.function <- match.arg(arg = combat.function)

  args.combat <- c('mod', 'par.prior', 'prior.plots', 'ref.batch', 'BBPARAM')
  args.combat_seq <- c('group', 'covar_mod', 'full_mod', 'shrink', 'shrink.disp',
                       'gene.subset.n')
  varargs <- list(...)

  layers <- layers %||% "data"
  scale.layer <- scale.layer %||% "scale.data"

  groups <- groups %||% Seurat:::CreateIntegrationGroups(object = object,
                                                         layers = layers,
                                                         scale.layer = scale.layer)
  groups.name <- groups.name %||% colnames(groups)[1]
  groups.name <- intersect(colnames(groups), groups.name)
  if (! length(x = groups.name)) {
    abort(message = "'groups.name' not in 'groups' data frame")
  }
  if (length(x = groups.name) > 1) {
    groups.name <- groups.name
    warning(paste("more 'groups.name' that expected. Using the first one",
                  sQuote(x = groups.name)), call. = FALSE, immediate. = TRUE)
  }

  if (use.scaled) {
    if (! scale.layer %in% Layers(object = object, search = scale.layer)) {
      abort(message = paste(sQuote(x = scale.layer), "not in object layers"))
    }
    layer.use <- scale.layer
  } else {
    layer.use <- unique(sub("\\..*", "", layers))
    if(length(layer.use) != 1) {
      abort(message="cannot find a consensus layer")
    }
    msg <- sprintf(fmt = "Layers used: %s\n",
                   paste(sQuote(Layers(object, search = layer.use)), collapse = ", "))
    message(msg[verbose], appendLF = FALSE)
    object <- JoinLayers(object = object, layers = layer.use)
  }
  data <- GetAssayData(object = object, layer = layer.use)
  data <- data[features %||% rownames(data), ]
  message(sprintf("Using %d features\n", nrow(data))[verbose], appendLF = FALSE)

  args <- list(batch = groups[colnames(data), groups.name, drop = TRUE])
  .ComBat <- ComBat
  if(combat.function == "combat") {
    args <- c(args, list(dat = data, mean.only = FALSE),
              varargs[intersect(names(varargs), args.combat)])
  } else {
    args <- c(args, list(counts = as.matrix(data)),
              varargs[intersect(names(varargs), args.combat_seq)])
    .ComBat <- ComBat_seq
  }
  corrected.mat <- do.call(what = .ComBat, args = args)
  rownames(corrected.mat) <- rownames(data)
  colnames(corrected.mat) <- colnames(data)
  output.list <- list()
  output.list[[reconstructed.assay]] <- CreateAssayObject(
    data = as(corrected.mat, "sparseMatrix"),
    key = key.assay)
  if(use.scaled) {
    output.list[[reconstructed.assay]] <- SetAssayData(
      object = output.list[[reconstructed.assay]],
      new.data = corrected.mat,layer = "scale.data")
  }
  return(output.list)

}

attr(x = CombatIntegration, which = 'Seurat.method') <- 'integration'
