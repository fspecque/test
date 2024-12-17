#' @include metrics_pca.R
#' @include utils.R
NULL

#' Score an embedding or a count matrix with the average silhouette width
#'
#' @description
#' Compute scores based on the average silhouette width (ASW) metric.
#'
#' \code{ScoreASW}: First, cell-to-cell distances are computed on the provided
#' embedding or layer matrix. Then, the silhouette score is calculated to
#' estimate the quality of the partition according to the variable with cell
#' type labels. Hence, this score measures to what extent cells with identical
#' label cluster together.
#'
#' \code{ScoreASWBatch}: Similar, but with the batch variable. This score
#' provides an estimation of batch mixing.
#'
#' @inheritParams integration-method
#' @param object A Seurat object
#' @param cell.var The name of the column with cell type label variable
#' (must be in the object metadata).
#' Ignored by \code{ScoreASWBatch} when \code{per.cell.var = FALSE}
#' @param what the name of the dimension reduction or layer to score. Must be
#' in the Seurat object or obtainable after a
#' \code{\link[SeuratObject:JoinLayers]{JoinLayers()}} call.
#' @param assay name of the assay to use. The output of
#' \code{\link[SeuratObject:DefaultAssay]{DefaultAssay()}} is used by default
#' @param metric name of the distance metric to use. One of 'euclidean',
#' 'cosine', 'angular', 'manhattan', 'hamming'. See \strong{Note} for details.
#' @param dist.package name of the package to compute distances with. One of
#' 'distances', 'Rfast' ,'parallelDist', 'stats'. The latter is always available,
#' the others must be installed beforehand. They are ordered from fastest to
#' slowest. When the requested package is not installed, the fastest amongst the
#' available ones is picked. See \strong{Note} for details.
#' @param ... additional parameters to pass to the distance computation
#' functions
#'
#' @return \code{ScoreASW} and \code{ScoreASWBatch}: a single float between 0
#' and 1, corresponding to the scaled average silhouette score.
#'
#' \code{AddScoreASW} and \code{AddScoreASWBatch}:  the updated Seurat
#' \code{object} with the ASW score(s) set for the integration.
#'
#' @importFrom SeuratObject DefaultAssay Reductions Embeddings Layers LayerData JoinLayers GetAssayData
#' @importFrom Matrix t
#' @importFrom rlang is_installed
#' @importFrom cluster silhouette
#'
#' @export
#' @details
#' \code{ScoreASW}: Given a matrix (reduction dimension or layer), the
#' cell-to-cell distance matrix \eqn{D} is computed. Then, the silhouette width
#' \eqn{s(i)} is calculated for each cell \eqn{i} with a label \eqn{c \in L}
#' (\eqn{L} is the set of possible cell type labels). Then, the mean of all
#' \eqn{s(i)} is computed (i.e. the ASW) and scaled between 0 and 1:
#' \deqn{\displaystyle ASW = \frac{1}{\left| L \right|} \times \sum_{i}{s(i)} \\[10pt]
#' score = \frac{ASW + 1}{2}}
#'
#' @note Those scores are an adaptation of the (cell-type) ASW and the batch ASW
#' as described in Luecken M.D. \emph{et al.}, 2022.
#'
#' Hamming distance is only supported by the \pkg{parallelDist} package, while
#' \pkg{distances} can only compute euclidean and related distances (cosine and
#' angular). Angular distances are actually refereed to as 'cosine' in
#' \code{\link[Seurat:FindNeighbors]{FindNeighbors()}} (\code{annoy.metric}),
#' hence called 'angular' here. Actual cosine dissimilarity-derived distances
#' are returned when \code{metric = 'cosine'}. Internally, angular distances are
#' computed with euclidean distances on \eqn{L^2} norm. cosine distances are
#' further transformed with :
#' \deqn{\displaystyle D_{cosine} = \frac{D_{angular}^2}{2}}
#'
#' @references Luecken, M. D., Büttner, M., Chaichoompu, K., Danese, A.,
#' Interlandi, M., Mueller, M. F., Strobl, D. C., Zappia, L., Dugas, M.,
#' Colomé-Tatché, M. & Theis, F. J. Benchmarking atlas-level data integration in
#' single-cell genomics. Nat Methods 19, 41–50 (2021).
#' \href{https://doi.org/10.1038/s41592-021-01336-8}{DOI}
#'
#' @seealso \code{\link[Seurat]{FindNeighbors}}, \code{\link[cluster]{silhouette}}
#' to know more about the silhouette metric, \code{\link{GetNeighborsPerBatch}}
#' and \code{\link{GetPropInterBatch}}
#'
#' @rdname score-asw
#' @name score-asw

ScoreASW <- function(object, cell.var,  what, assay = NULL,
                     metric = c('euclidean', 'cosine', 'angular', 'manhattan', 'hamming'),
                     dist.package = c('distances', 'Rfast' ,'parallelDist', 'stats'),
                     verbose = TRUE, ...) {
  metric <- tolower(metric)
  metric <- match.arg(metric)
  metric_ <- ifelse(metric %in% c('cosine', 'angular'), 'euclidean', metric)
  dist.package <- match.arg(dist.package)
  assay <- assay %||% DefaultAssay(object)

  df.mtdt <- object[[]]
  cell.var.in <- cell.var %in% colnames(df.mtdt)
  msg <- "are absent from colnames of metadata"
  if (sum(cell.var.in) == 0) {
    abort(paste("All the provided cell.var", msg))
  }
  if (sum(cell.var.in) < length(cell.var)) {
    warning(sprintf("%d out of %d cell.var %s (%s). Ignoring them",
                    sum(!cell.var.in), length(cell.var), msg,
                    paste(sQuote(cell.var[!cell.var.in]), collapse = ', ')),
            call. = F, immediate. = T)
    cell.var <- cell.var[cell.var.in]
  }
  cell.var <- lapply(asplit(df.mtdt[, cell.var, drop = FALSE], MARGIN = 2),
                     function(vec) as.integer(as.factor(vec)))

  if (what %in% Reductions(object)) {
    mat <- Embeddings(object, reduction = what)
    message(sprintf('found %s in Seurat object\'s reductions\n', sQuote(what))[verbose],
            appendLF = FALSE)
  } else if (what %in% Layers(object)) {
    mat <- t(LayerData(object, layer = what, assay = assay))
    message(sprintf('found %s in Seurat object\'s layers\n', sQuote(what))[verbose],
            appendLF = FALSE)
  } else {
    l <- Layers(object, search = what, assay = assay)
    if (l %iff% TRUE %||% FALSE) {
      if (inherits(sub.object[[assay]], "StdAssay")) {
        mat <- t(GetAssayData(JoinLayers(object, assay = assay, layers = what),
                              assay = assay, layer = what))
      } else {
        mat <- t(GetAssayData(object, assay = assay, layer = what))
      }

      message(sprintf('found %s in Seurat object\'s layers after join\n', sQuote(what))[verbose],
              appendLF = FALSE)
    } else {
      abort(paste(sQuote(what), 'does not seem to be a reduction nor a layer'))
    }
  }

  if (dist.package != 'stats' && ! is_installed(dist.package)) {
    possible.pkg <- c('distances', 'Rfast' ,'parallelDist', 'stats')
    installed.pkg <- sapply(possible.pkg, is_installed)
    warning(sQuote(dist.package), ' not installed, using ',
            sQuote(possible.pkg[installed.pkg][[1]]), ' instead.',
            call. = F, immediate. = T)
    dist.package <- possible.pkg[installed.pkg][[1]]
  }
  if (metric_ == 'hamming' && dist.package != 'parallelDist') {
    abort('hamming distance is only supported by parallelDist')
  } else if (dist.package == "distances" && metric  == 'manhattan') {
    abort('distances package does not support manhattan distances')
  }

  if (metric %in% c('cosine', 'angular')) { # L2 norm (M / rowwise(sqrt(sum(M^2))))
    mat <- NormaliseL2(mat = mat, MARGIN = 1L)
    # euclidean distance on mat -> ~ 'angular', i.e. 'cosine' of FindNeighbors (AnnoyAngular)
  }

  mat <- as.matrix(mat)
  message(sprintf('Computing %s distances using %s...', metric,
                  sQuote(dist.package))[verbose], appendLF = FALSE)
  dist.mat <- switch (dist.package,
                      distances = distances::distance_matrix(distances::distances(mat), ...),
                      Rfast = Rfast::Dist(mat, method = metric_, ...),
                      parallelDist = parallelDist::parallelDist(mat, method = metric_, ...),
                      stats = stats::dist(mat, method = metric_, ...)
  )
  if (metric == 'cosine') {
    dist.mat <- dist.mat**2 / 2
  }
  message('done\n'[verbose], appendLF = FALSE)

  sils <- sapply(cell.var, function(vec)
    silhouette(x = vec, dist = dist.mat)[,'sil_width'],
                 simplify = FALSE, USE.NAMES = TRUE)
  scores <- sapply(sils, function(sil) (mean(sil) + 1) / 2,
                   simplify = TRUE, USE.NAMES = TRUE)
  return(scores)
}

#' @param integration name of the integration to score
#' @export
#' @rdname score-asw
AddScoreASW <- function(object, integration,
                        cell.var,  what, assay = NULL,
                        metric = c('euclidean', 'cosine', 'angular', 'manhattan', 'hamming'),
                        dist.package = c('distances', 'Rfast' ,'parallelDist', 'stats'),
                        verbose = TRUE, ...) {
  scores <- ScoreASW(object, cell.var = cell.var, what = what, assay = assay,
                     metric = metric, dist.package = dist.package,
                     verbose = verbose, ...)

  score.names <- paste("ASW", names(scores), sep = '_')
  object <- check_misc(object)
  for (i in 1:length(scores)) {
    object <- SetMiscScore(object, integration = integration,
                           score.name = score.names[i],
                           score.value = scores[[i]])
  }
  return(object)
}

#' @inheritParams ScoreASW
#' @param batch.var The name of the column with batch variable.
#' (must be in the object metadata). Required by \code{ScoreASWBatch}.
#' @param per.cell.var whether to compute silhouette coefficients with the batch
#' variable for each cell-type separately (default behaviour). Setting to
#' \code{FALSE} causes the silhouette coefficients to be computed on the whole
#' data directly.
#'
#' @importFrom SeuratObject DefaultAssay Reductions Embeddings Layers LayerData JoinLayers GetAssayData
#' @importFrom Matrix t
#' @importFrom rlang is_installed :=
#' @importFrom cluster silhouette
#'
#' @export
#' @details
#' \code{ScoreAWSBatch}: The default parameters (\code{per.cell.var = TRUE})
#' correspond to the original score from Luecken M.D. \emph{et al.}, 2022. It is
#' computed as follow: for each cell type label \eqn{c} among the set of all
#' labels \eqn{L}, the cell-to-cell matrix distance \eqn{D_c} is computed for
#' cells \eqn{i} with label \eqn{c}. Each cell's silhouette width \eqn{s(i)} is
#' calculated according to the batch variable and transformed to be as close to
#' 1 as its absolute value is close to 0. An average silhouette width
#' \eqn{ASW_c} is then computed per cell type label \eqn{c \in L} and the mean
#' of those correspond to the final score:
#' \deqn{\displaystyle ASW_c = \frac{1}{\left| c \right|} \times \sum_{i \in c}{1 - \left| s(i) \right|} \\[10pt]
#' score = \frac{1}{\left| L \right|} \times \sum_{c \in L}{ASW_c}}
#'
#' When \code{per.cell.var = FALSE}, \eqn{ASW} is computed for all cells at once
#' (just like for \code{ScoreASW} but on the batch variable), then scaled
#' and averaged similarly:
#' \deqn{\displaystyle score = ASW = \frac{1}{N} \times \sum_{i=1}^{N}{1 - \left| s(i) \right|}}
#' with \eqn{N} being the total number of cells
#' @rdname score-asw
#' @name score-asw
ScoreASWBatch <- function(object, batch.var = NULL, cell.var = NULL,  what,
                          per.cell.var = TRUE, assay = NULL,
                          metric = c('euclidean', 'cosine', 'angular', 'manhattan', 'hamming'),
                          dist.package = c('distances', 'Rfast' ,'parallelDist', 'stats'),
                          verbose = TRUE, ...) {
  metric <- tolower(metric)
  metric <- match.arg(metric)
  metric_ <- ifelse(metric %in% c('cosine', 'angular'), 'euclidean', metric)
  dist.package <- match.arg(dist.package)
  assay <- assay %||% DefaultAssay(object)

  .prep_MetaDataBatch(object = object, batch.var = batch.var,
                      assay = assay, layer = 'data')
  if (per.cell.var) {
    cell.var.in <- cell.var %in% colnames(df.mtdt)
    msg <- "are absent from colnames of metadata"
    if (sum(cell.var.in) == 0) {
      abort(paste("All the provided cell.var", msg))
    }
    if (sum(cell.var.in) < length(cell.var)) {
      warning(sprintf("%d out of %d cell.var %s (%s). Ignoring them",
                      sum(!cell.var.in), length(cell.var), msg,
                      paste(sQuote(cell.var[!cell.var.in]), collapse = ', ')),
              call. = F, immediate. = T)
      cell.var <- cell.var[cell.var.in]
    }
  } else {
    cell.var <- 'identity'
    df.mtdt$identity <- 1
  }

  if (what %in% Reductions(object)) {
    mat <- Embeddings(object, reduction = what)
    message(sprintf('found %s in Seurat object\'s reductions\n', sQuote(what))[verbose],
            appendLF = FALSE)
  } else if (what %in% Layers(object)) {
    mat <- t(LayerData(object, layer = what, assay = assay))
    message(sprintf('found %s in Seurat object\'s layers\n', sQuote(what))[verbose],
            appendLF = FALSE)
  } else {
    l <- Layers(object, search = what, assay = assay)
    if (l %iff% TRUE %||% FALSE) {
      mat <- t(GetAssayData(JoinLayers(object, assay = assay, layers = what),
                            assay = assay, layer = what))
      message(sprintf('found %s in Seurat object\'s layers after join\n', sQuote(what))[verbose],
              appendLF = FALSE)
    } else {
      abort(paste(sQuote(what), 'does not seem to be a reduction nor a layer'))
    }
  }

  if (dist.package != 'stats' && ! is_installed(dist.package)) {
    possible.pkg <- c('distances', 'Rfast' ,'parallelDist', 'stats')
    installed.pkg <- sapply(possible.pkg, is_installed)
    warning(sQuote(dist.package), ' not installed, using ',
            sQuote(possible.pkg[installed.pkg][[1]]), ' instead.',
            call. = F, immediate. = T)
    dist.package <- possible.pkg[installed.pkg][[1]]
  }
  if (metric_ == 'hamming' && dist.package != 'parallelDist') {
    abort('hamming distance is only supported by parallelDist')
  } else if (dist.package == "distances" && metric  == 'manhattan') {
    abort('distances package does not support manhattan distances')
  }

  if (metric %in% c('cosine', 'angular')) { # L2 norm (M / rowwise(sqrt(sum(M^2))))
    mat <- NormaliseL2(mat = mat, MARGIN = 1L)
    # euclidean distance on mat -> ~ 'angular', i.e. 'cosine' of FindNeighbors (AnnoyAngular)
  }

  scores <- sapply(cell.var, function(vec) {
    i <- 0
    prog.bar <- NULL
    if (verbose) {
      n <- length(unique(df.mtdt[, cell.var, drop = TRUE]))
      if (n > 1) {
        prog.bar <- utils::txtProgressBar(i, n, width = 30, style = 3)
      }
    }
    sils <- df.mtdt %>%
      mutate(!!sym(batch.var) := as.integer(as.factor(!!sym(batch.var)))) %>%
      group_by(!!sym(vec)) %>% group_map(
        function(df.mtdt.sub, y) {
          i <<- i + 1
          prog.bar %iff% setTxtProgressBar(prog.bar, i)
          batches <- df.mtdt.sub[, batch.var, drop = TRUE]
          if (length(unique(batches)) %in% c(1, nrow(df.mtdt.sub))) return(NULL)

          cells <- df.mtdt.sub[, idcol, drop = TRUE]
          sub.mat <- as.matrix(mat[cells, , drop = FALSE])
          dist.mat <- switch (
            dist.package,
            distances = distances::distance_matrix(distances::distances(sub.mat)),
            Rfast = Rfast::Dist(sub.mat, method = metric_),
            parallelDist = parallelDist::parallelDist(sub.mat, method = metric_),
            stats = stats::dist(sub.mat, method = metric_)
          )
          if (metric == 'cosine') {
            dist.mat <- dist.mat**2 / 2
          }

          sil <- silhouette(x = df.mtdt.sub[, batch.var, drop = TRUE],
                            dist = dist.mat)[,'sil_width']
          return(mean(abs(sil)))
        }
      )
    prog.bar %iff% close(prog.bar)

    sils <- unlist(sils, use.names = FALSE) %||% NA
    return(mean(1 - sils))
  }, simplify = "numeric", USE.NAMES = TRUE)
  return(scores)
}

#' @param integration name of the integration to score
#' @export
#' @rdname score-asw
AddScoreASWBatch <- function(object, integration,
                             batch.var = NULL, cell.var = NULL,  what,
                             per.cell.var = TRUE, assay = NULL,
                             metric = c('euclidean', 'cosine', 'angular', 'manhattan', 'hamming'),
                             dist.package = c('distances', 'Rfast' ,'parallelDist', 'stats'),
                             verbose = TRUE, ...) {
  scores <- ScoreASWBatch(object, batch.var = batch.var, cell.var = cell.var,
                          what = what, per.cell.var = per.cell.var,
                          assay = assay, metric = metric,
                          dist.package = dist.package, verbose = verbose, ...)

  score.names <- paste("ASW.batch", names(scores), sep = '_')
  object <- check_misc(object)
  for (i in 1:length(scores)) {
    object <- SetMiscScore(object, integration = integration,
                           score.name = score.names[i],
                           score.value = scores[[i]])
  }
  return(object)
}
