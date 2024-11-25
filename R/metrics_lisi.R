#' @include utils.R
#' @include expand_knn_graph.R
NULL

#' Score a dimensionality reduction embedding or knn graph using the Local
#' Inverse Simpson Index
#'
#' @description
#' Compute the Local Inverse Simpson's Index (LISI) to estimate batch mixing or
#' cell type mixing (iLISI and cLISI respectively according to Luecken M.D.
#' \emph{et al.}, 2022).
#'
#' \code{AddScoreLISI} returns an updated Seurat object, while \code{ScoreLISI}
#' outputs the raw LISI scores for each cell
#'
#' @param object A Seurat object
#' @param batch.var The name of the batch variable (must be in the object metadata).
#' Can be omitted if \code{cell.var} is not \code{NULL}
#' @param cell.var The name of the  cell variable (must be in the object metadata).
#' Can be omitted if \code{batch.var} is not \code{NULL}
#' @param reduction The name of the reduction to score. Arguments
#' \code{reduction} and \code{graph.name} are mutually exclusive
#' @param dims The dimensions of \code{reduction} to consider. All dimensions are
#' used by default. Ignored when scoring a graph
#' @param graph.name The name of the knn graph to score. Arguments
#' \code{reduction} and \code{graph.name} are mutually exclusive
#' @param graph.type one of 'distances' or 'connectivities' (not supported yet).
#' Ignored when scoring a cell embedding
#' @param do.symmetrize whether to symmetrize the knn graphs. Set to\code{FALSE}
#' to disable (not recommended, especially when scoring a knn graph directly.
#' See \strong{Details})
#' @param save.graph whether to save the graph used to compute the LISI score(s)
#' in the Seurat object
#' @param return.graph whether to return the graph used to compute the LISI score(s)
#' @param new.graph name of the graph used to compute the LISI score(s).
#' When \code{new.graph = NULL} (default), a name is constructed depending on
#' input and arguments. Ignored when \code{save.graph = FALSE}
#' @param perplexity third of number of each cell's neighbours. When the value of
#' perplexity and the number of neighbours (*3) are discrepant, the graph is
#' adapted. Multiple scores with varying values of perplexity are not comparable,
#' hence it is recommended to use the same value for each integration to score.
#' @param tol Stop the computation of the local Simpson's index when it converges
#' to this tolerance.
#' @param do.scale whether to scale the output LISI values between 0 and 1.
#' @param largest.component.only whether to compute the scores on the largest
#' component or all sufficiently large components (default)
#' @param assay the name of the assay to reference in the output Graph object (
#' when \code{save.graph = TRUE})
#' @param verbose whether to print progress messages
#' @param ... Additional parameters to pass to other methods
#' (see \strong{Details}).
#'
#' @return \code{ScoreLISI}: a data frame with each cell's raw LISI score, or
#' a list containing the aforementioned data frame and the graph used to compute
#' it (\code{return.graph = TRUE}).
#'
#' \code{AddScoreLISI}: the updated Seurat \code{object}, with cell-wise LISI
#' scores in the meta data (identical to \code{ScoreLISI}'s output), global
#' scores in misc and a new Graph object when \code{save.graph = TRUE}.
#'
#' @importFrom rlang maybe_missing
#' @importFrom SeuratObject AddMetaData Misc<-
#' @importFrom dplyr %>% rename all_of summarize across n_distinct cur_column
#'
#' @export
#' @details
#' When scoring a reduction, a knn graph with enough neighbours per cell is
#' computed. If \code{do.symmetrize = TRUE}, the graph is symmetrized and the k
#' best neighbours are kept.
#'
#' When scoring a knn graph, the graph is expanded with Dijkstra's algorithm to
#' reach enough neighbours per cell. If \code{do.symmetrize = TRUE}, the graph
#' is symmetrized beforehand. Note that when \code{do.symmetrize = FALSE},
#' Dijkstra's algorithm is applied on the asymmetric distance matrix, and
#' components are computed. But each component is then symmetrized and Dijkstra's
#' algorithm is computed on each of them. Indeed, there is no guarantee that the
#' cells' number of neighbours is preserved after decomposing the directed graph.
#' Hence, when \code{do.symmetrize = FALSE} and a graph is input, the graph is
#' considered as directed only to find components.
#'
#' In either case, it is recommended to keep \code{do.symmetrize = TRUE}.
#'
#' For possible additional parameters, see \code{\link[Seurat]{FindNeighbors}}
#' (when inputting a reduction) or \code{\link{ExpandNeighbors}} (when inputting
#' a knn graph)
#'
#' @note This score is an adaptation of the LISI score as described in Korsunsky
#' I. \emph{et al.}, 2019 and also used in Luecken M.D. \emph{et al.}, 2022.
#'
#'
#' @references Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F.,
#' Wei, K., Baglaenko, Y., Brenner, M., Loh, P. & Raychaudhuri, S. Fast,
#' sensitive and accurate integration of single-cell data with Harmony. Nat
#' Methods 16, 1289–1296 (2019).
#' \href{https://doi.org/10.1038/s41592-019-0619-0}{DOI}
#' @references Luecken, M. D., Büttner, M., Chaichoompu, K., Danese, A.,
#' Interlandi, M., Mueller, M. F., Strobl, D. C., Zappia, L., Dugas, M.,
#' Colomé-Tatché, M. & Theis, F. J. Benchmarking atlas-level data integration in
#' single-cell genomics. Nat Methods 19, 41–50 (2021).
#' \href{https://doi.org/10.1038/s41592-021-01336-8}{DOI}
#'
#' @seealso \code{\link[Seurat]{FindNeighbors}}, \code{\link{ExpandNeighbors}}
#' and \code{\link{CutKnn}}
#' @rdname score-lisi

AddScoreLISI <- function(object, integration,
                         batch.var = NULL, cell.var = NULL,
                         # DimReduc
                         reduction, dims = NULL,
                         # Graph | Neighbor
                         graph.name,
                         graph.type = c("distances", "connectivities"),
                         do.symmetrize = TRUE,

                         save.graph = TRUE,
                         new.graph = NULL,
                         perplexity = 30, tol = 1e-5, do.scale = TRUE,
                         largest.component.only = FALSE,
                         assay = NULL,
                         verbose = TRUE, ...) {
  reduction <- maybe_missing(reduction, default = NULL)
  reduction <- if(isFALSE(reduction)) NULL else reduction
  graph.name <- maybe_missing(graph.name, default = NULL)
  graph.name <- if(isFALSE(graph.name)) NULL else graph.name
  object.name <- reduction %||% graph.name

  # i <- c(batch.var %iff% T %||% F, cell.var %iff% T %||% F)
  v <- c(batch.var, cell.var)
  # new.names <- paste0(c("LISIbatch_", "LISIcell_"), v, "_", object.name)[i]
  nb <- length(batch.var)
  nc <- length(cell.var)
  new.names <- paste0(rep(c("iLISI_", "cLISI_"), times = c(nb, nc)),
                      v)[seq_len(nb + nc)]

  lisi.out <- ScoreLISI(object = object, batch.var = batch.var,
                        cell.var = cell.var, reduction = reduction,
                        dims = dims, graph.name = graph.name,
                        graph.type = graph.type, do.symmetrize = do.symmetrize,
                        return.graph = TRUE, perplexity = perplexity,
                        tol = tol,
                        largest.component.only = largest.component.only,
                        assay = assay, verbose = verbose, ...)

  new.names <- setNames(v, new.names)
  new.names <- new.names[v %in% colnames(lisi.out$lisi)]
  v <- unname(new.names)
  n <- object[[]] %>% summarize(across({{ v }}, n_distinct))
  lisi.score <- lisi.out$lisi %>%
    mutate(across({{ v }}, ~ numeric_lisi(.x, N = n[[cur_column()]])))

  object <- AddMetaData(object = object,
                        metadata = lisi.score %>% rename(all_of(new.names)))
  if (save.graph) {
    new.graph <- new.graph %||%
      paste0(object.name, '_', 'symmetric_'[do.symmetrize],
             sprintf('dijkstra_%dk_', as.integer(3 * perplexity))[graph.name %iff% 1 %||% 0],
            'LISI_perp.', perplexity)
    object[[new.graph]] <- lisi.out$graph
    message(sprintf('New graph used for LISI saved as %s\n',
                    sQuote(new.graph))[verbose], appendLF = F)
  }

  n <- n %>% rename(all_of(new.names)) %>% unlist()
  lisi.score <- lisi.score %>%
    summarize(across({{ v }}, ~ median(.x))) %>%
    rename(all_of(new.names))
  object <- check_misc(object)
  for (idx in 1:length(lisi.score)) {
    object <- SetMiscScore(object, integration = integration,
                           score.name = names(lisi.score)[idx],
                           score.value = lisi.score[[idx]],
                           class = "numeric_lisi")
    # N(object@misc$si_scores[[names(lisi.score)[idx]]]) <- n[[idx]]
  }

  # if (do.scale) {
  #   lisi.score[batch.var] <- (lisi.score[batch.var] - 1) / (n[batch.var] - 1)
  #   lisi.score[cell.var] <- (n[cell.var] - lisi.score[cell.var]) / (n[cell.var] - 1)
  # }
  # slot <- sprintf("LISI_perp.%s_%s_%s", as.integer(perplexity),
  #                 reduction %||% graph.name, c("unscaled", "scaled")[do.scale + 1])
  # Misc(object = object, slot = slot) <- lisi.score %>% rename(all_of(new.names))

  object
}

#' @importFrom rlang maybe_missing
#' @importFrom SeuratObject Embeddings Graphs Neighbors
#' @importFrom dplyr %>% mutate across
#' @export
#' @rdname score-lisi

ScoreLISI <- function(object, batch.var = NULL, cell.var = NULL,
                      # DimReduc
                      reduction, dims = NULL,
                      # Graph | Neighbor
                      graph.name,
                      graph.type = c("distances", "connectivities"),
                      do.symmetrize = TRUE,

                      return.graph = FALSE,
                      perplexity = 30, tol = 1e-5,
                      largest.component.only = FALSE,
                      assay = NULL,
                      verbose = TRUE, ...) {
  reduction <- maybe_missing(reduction, default = NULL)
  reduction <- if(isFALSE(reduction)) NULL else reduction
  graph.name <- maybe_missing(graph.name, default = NULL)
  graph.name <- if(isFALSE(graph.name)) NULL else graph.name

  reduction %||% graph.name %||%
    abort('a graph or a reduction to score is needed but both are missing')
  reduction %iff% graph.name %iff%
    abort('you cannot score a reduction and a graph at the same time')

  graph.type <- tolower(graph.type)
  graph.type <- match.arg(graph.type)
  ##############################################################################
  if (graph.type == "connectivities") {
    abort("not implemented")
  }
  ##############################################################################
  return.graph <- return.graph %||% FALSE

  # return.graph <- list(...)[["return.graph"]] %||% FALSE
  has.vars <- c("batch" = !is.logical(batch.var %||% FALSE),
                "cell"  = !is.logical(cell.var  %||% FALSE))
  if (!any(has.vars)) {
    abort(message = sprintf("please specify at lease one of %s or %s",
                            sQuote("batch.var"), sQuote("cell.var")))
  }
  vars <- c("batch" = batch.var, "cell" = cell.var)
  mtdt <- object[[]]
  found.vars <- intersect(vars, colnames(mtdt))
  if (length(found.vars) == 0) {
    abort(message = sprintf("neither %s were found in the colnames of the metadata",
                            paste(sQuote(vars),
                                  collapse = " nor ")))
  }
  if (length(found.vars) < length(vars)) {
    warning(sprintf("%s was not found in the colnames of the metadata, skipping",
                    sQuote(setdiff(vars, found.vars))),
            call. = FALSE, immediate. = TRUE)
  }
  if (reduction %iff% TRUE %||% FALSE) {
    if (! reduction %in% Reductions(object)) {
      abort(message = sprintf("%s reduction not in object", sQuote(reduction)))
    }
    object_ <- Embeddings(object = object, reduction = reduction)
  } else { #graph
    if (! graph.name %in% c(Graphs(object), Neighbors(object))) {
      abort(message = sprintf("%s graph not in object", sQuote(graph.name)))
    }
    object_ <- suppressWarnings(
      Graphs(object, graph.name) %||% Neighbors(object, graph.name))
  }
  graph.objects <- .prep.lisi(object = object_,
                              dims = dims, graph.type = graph.type,
                              do.symmetrize = do.symmetrize,
                              perplexity = perplexity, tol = tol,
                              largest.component.only = largest.component.only,
                              assay = assay,
                              verbose = verbose, ...)


  lisi.installed <- requireNamespace("lisi", quietly = TRUE)
  mtdt <- mtdt %>% mutate(across({{ found.vars }}, ~ as.integer(as.factor(.x)) -
                                   as.integer(lisi.installed)))
  if (! lisi.installed) {
    FUN.LSI <- .compute.lsi
    if (pb$print.lisi.msg) {
      message("To benefit from the original and much faster implementation of ",
              "the Local Inverse Simpson’s Index (LISI), you can install the ",
              "lisi package (by the same team that developed harmony) using either:",
              "\n\nremotes::install_github('immunogenomics/lisi')",
              "\n# or",
              "\ndevtools::install_github('immunogenomics/lisi')",
              "\n\nThis message will be shown once per session")
      pb$print.lisi.msg <- FALSE
    }
    message("Computing LISI scores..."[verbose], appendLF = F)
  } else {
    FUN.LSI <- lisi::compute_simpson_index
    message("Computing LISI scores with lisi package..."[verbose], appendLF = F)
  }

  simpson.indexes <- lapply(graph.objects[["compos"]], function(graph.object) {
    nn.idx  <- t(slot(graph.object, 'nn.idx')[,-1]) - as.integer(lisi.installed)
    nn.dist <- t(slot(graph.object, 'nn.dist')[,-1])
    cell.names <- slot(graph.object, 'cell.names')
    si <- lapply(asplit(mtdt[, found.vars, drop = FALSE], 2), function(v) {
      FUN.LSI(nn.dist, nn.idx, v[cell.names], perplexity = perplexity, tol = tol,
              n_batches = length(unique(v)))
    })

    si <- do.call(cbind, si)
    rownames(si) <- cell.names
    colnames(si) <- found.vars
    si
  })
  message("done.\n"[verbose], appendLF = F)
  simpson.indexes <- 1 / as.data.frame(do.call(rbind, simpson.indexes)) # inverse SI

  if (return.graph) {
    simpson.indexes <- list(graph = graph.objects[["full"]], lisi = simpson.indexes)
  }
  return(simpson.indexes)
}

.prep.lisi <- function(object, dims = NULL,
                       graph.type = c("distances", "connectivities"),
                       do.symmetrize = TRUE, perplexity = 30, tol = 1e-5,
                       largest.component.only = FALSE,
                       assay = NULL,
                       verbose = TRUE, ...) {
  UseMethod(generic = '.prep.lisi', object = object)
}

#' @importFrom Seurat FindNeighbors
.prep.lisi.DimReduc <- function(object, dims = NULL, do.symmetrize = TRUE,
                                perplexity = 30, tol = 1e-5,
                                largest.component.only = FALSE,
                                graph.type = c("distances", "connectivities"),
                                assay = NULL,
                                verbose = TRUE, ...) {
  dims <- dims %||% 1:ncol(object)
  k.lisi <- as.integer(3*perplexity)
  dot.args <- list(...)
  if (! all(dims %in% 1:ncol(object))) {
    abort(message = "some dims are out of range")
  }
  object <- object[, dims, drop = FALSE]  # output matrix if DimReduc
  fixed.args <- list(object = object, query = NULL, distance.matrix = FALSE,
                     k.param = k.lisi,
                     return.neighbor = TRUE, compute.SNN = FALSE,
                     verbose = verbose)
  default.args <- formals('FindNeighbors.default', envir = asNamespace('Seurat'))
  names.args <- setdiff(names(default.args), c(names(fixed.args), '...'))
  default.args <- Map(`%||%`, x = dot.args[names.args],
                      y = default.args[names.args])
  default.args <- setNames(default.args, names.args)
  args <- c(fixed.args, default.args)

  graph.object <- do.call(FindNeighbors, args)

  return(.prep.lisi(object = graph.object, do.symmetrize = do.symmetrize,
                    graph.type = "dist", perplexity = perplexity, tol = tol,
                    largest.component.only = largest.component.only,
                    from.dimred = TRUE, assay = assay, verbose = verbose, ...))
}

.prep.lisi.matrix <- .prep.lisi.DimReduc

#' @importFrom SeuratObject as.Graph
#' @importFrom igraph components graph_from_adjacency_matrix
#' @importFrom Matrix drop0 diag<-
.prep.lisi.Graph <- function(object, do.symmetrize = TRUE,
                             graph.type = c("distances", "connectivities"),
                             perplexity = 30, tol = 1e-5,
                             largest.component.only = FALSE,
                             from.dimred = FALSE,
                             dims = NULL,
                             assay = NULL,
                             verbose = TRUE, ...) {
  graph.type <- tolower(graph.type)
  graph.type <- match.arg(graph.type)
  k.lisi <- as.integer(3*perplexity)
  dot.args <- list(...)
  if (do.symmetrize) {
    object <- SymmetrizeKnn(object, use.max = TRUE)
  }
  if (! from.dimred && get.k(object, FUN = 'min') < k.lisi) {
    fixed.args <- list(object = object, graph.type = graph.type,
                       k.target = k.lisi,
                       do.symmetrize = do.symmetrize,
                       verbose = verbose)
    default.args <- formals('expand_neighbours_dijkstra',
                            envir = asNamespace(.packageName))
    names.args <- setdiff(names(default.args), c(names(fixed.args), '...'))
    default.args <- Map(`%||%`, x = dot.args[names.args],
                        y = default.args[names.args])
    default.args <- setNames(default.args, names.args)
    args <- c(fixed.args, default.args)
    args$which.dijkstra <- match.arg(eval(args$which.dijkstra),
                                     eval(default.args$which.dijkstra))
    object <- do.call(expand_neighbours_dijkstra, args)
  }

  k.graph <- get.k(object, FUN = 'max')
  if (k.graph < k.lisi) {
    abort(sprintf('cannot gather enough neighbours (%d at beast, %d required)
Try to decrease `perplexity` to %d or lower', k.graph, k.lisi, k.graph %/% 3))
  }
  graph.object <- .CutKnn(object = object, k.max = k.lisi,
                          assay = assay, verbose = verbose)

  if (!inherits(graph.object, "Matrix")) {
    graph.object <- as.Graph(graph.object)
  }
  diag(graph.object) <- 1
  mode_graph <- paste0('un'[do.symmetrize], 'directed')
  mode_compo <- c('strong', 'weak')[[do.symmetrize + 1]]
  graph.compos <- components(
    graph_from_adjacency_matrix(graph.object > 0, mode = mode_graph,
                                weighted = NULL,diag = TRUE),
    mode = mode_compo
  )
  k.thresh <- c(k.lisi, max(graph.compos$csize))[[largest.component.only + 1]]
  compos.keep <- which(graph.compos$csize >= k.thresh)
  # }
  idx <- which(graph.compos$membership %in% compos.keep)
  cells.keep <- graph.compos$membership[idx]
  if (! do.symmetrize & ! from.dimred) {
    graph.objects <- lapply(split(names(cells.keep), cells.keep), function(ij)
      as.Graph(drop0(graph.object[ij, ij, drop=F])))
    graph.objects <- lapply(graph.objects, function(x) {
      .prep.lisi.Graph(x, do.symmetrize = TRUE, graph.type = graph.type,
                       perplexity = perplexity, tol = tol,
                       largest.component.only = T, from.dimred = FALSE,
                       assay = assay, verbose = FALSE, ...)$compos[[1]]
    })
  } else {
    graph.objects <- lapply(split(names(cells.keep), cells.keep), function(ij) {
      g <- graph.object[ij, ij, drop=F]
      diag(g) <- -Inf
      as.Neighbor(as.Graph(g))
    })
  }
  diag(graph.object) <- 0
  return(list(full = as.Graph(drop0(graph.object)), compos = graph.objects))
}

.prep.lisi.Neighbor <- function(object, do.symmetrize = TRUE,
                                graph.type = c("distances", "connectivities"),
                                perplexity = 30, tol = 1e-5,
                                largest.component.only = FALSE,
                                from.dimred = FALSE,
                                assay = NULL,
                                verbose = TRUE, ...) {
  object <- as.Graph(object)
  return(.prep.lisi(object = object, do.symmetrize = do.symmetrize,
                    graph.type = graph.type, perplexity = perplexity, tol = tol,
                    largest.component.only = largest.component.only,
                    from.dimred = from.dimred, assay = assay, verbose = verbose,
                    ...))
}


.compute.lsi <- function(nn.dist, nn.idx, batch.var, perplexity = 30, tol = 1e-5, ...) {
  ncells = ncol(nn.dist)
  nbatch = list(...)$n_batches %||% length(unique(batch.var))
  proba = rep(0, nrow(nn.dist))
  simpson = rep(0, ncells)
  logU = log(perplexity)

  for (i in 1:ncells) {
    beta <- 1
    betamin <- -Inf
    betamax <- Inf
    cell.nn.dist <- nn.dist[, i, drop=TRUE]
    e <- environment()
    list2env(.Hbeta(cell.nn.dist, beta, proba), e)
    Hdiff <- H - logU
    tries <- 0
    # Step 1: get neighbour probabilities
    while(abs(Hdiff) > tol && tries < 50) {
      if (Hdiff > 0){
        betamin <- beta
        if (is.infinite(betamax)) {beta <- beta * 2}
        else {beta <- (beta + betamax) / 2}
      } else{
        betamax <- beta
        if (is.infinite(betamin)) {beta <- beta /2}
        else {beta <- (beta + betamin) / 2}
      }

      list2env(.Hbeta(cell.nn.dist, beta, proba), e)
      Hdiff <- H - logU
      tries <- tries + 1
    }

    if (H == 0) {
      simpson[i] <- -1
    }

    # Step 2: compute Simpson's Index
    for (b in 1:nbatch) {
      q <- which(batch.var[nn.idx[,i]] == b)  # indices of cells belonging to batch (b)
      if (length(q) > 0) {
        sumP <- sum(proba[q])
        simpson[i] <- simpson[i] + sumP^2
      }
    }
  }
  return(simpson)
}

.Hbeta <- function(cell.nn.dist, beta, proba) {
  proba <- exp(- cell.nn.dist * beta)
  sumP <- sum(proba)
  if (sumP == 0) {
    H <- 0
    proba <- rep(0, length(cell.nn.dist))
  } else {
    H <- log(sumP) + beta * sum(cell.nn.dist * proba) / sumP
    proba <- proba / sumP
  }
  return(list("H" = H, "proba" = proba))
}
