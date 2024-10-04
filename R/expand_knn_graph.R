#' @include utils.R
NULL

#' Expand knn graph to increase the number of neighbours
#'
#' @description
#' Expand a knn graph to increase the number of nearest neighbours using
#' Dijkstra's algorithm, or a diffusion algorithm. Dijkstra's algorithm is used
#' to prepare for LISI score, while diffusion is suited for preparing to compute
#' the kBET score. Beware that diffusion is designed to work on connectivity
#' matrices and is not adequate for distance-based networks.
#'
#' @inheritParams SymmetrizeKnn
#' @param new.graph.name name of the return \code{Graph} to store in the
#' \code{Seurat object}.
#' @param k.target number of nearest neighbours to reach
#' @param do.symmetrize whether to make the input graph symmetric if necessary.
#' See \strong{Details} section for further explanations
#' @param algo One of "dijkstra" or "diffusion". "diffusion" is suited for
#' connectivity matrices only
#' @param which.dijkstra one of "igraph", "fast" or "slow". "auto" (default)
#' chooses for you. See \strong{Details} section
#' @param dijkstra.ncores number of cores to use for Dijkstra's algorithm.
#' Ignored when \code{which.dijkstra = "igraph"}
#' @param dijkstra.tol number of sequential iterations with identical best
#' neighbours found to consider that Dijkstra's algorithm should be stopped.
#' Ignored when \code{which.dijkstra = "igraph"}
#' @param diffusion.iter maximum number of iterations to reach \code{k.target}
#' @param verbose whether to print progress messages
#'
#' @return the Seurat object with a new \code{Graph} instance or a
#' \code{dgCMatrix} representing the \code{Graph} itself
#'
#' @details
#' The approximate nature of the nearest neighbour search algorithm used to
#' compute the knn graph makes the resulting adjacency matrix asymmetric.
#' It means that \eqn{cell_i} can be a nearest neighbour of \eqn{cell_j} while
#' the reverse is not true even though the distance between \eqn{cell_i} and
#' \eqn{cell_j} is lower than the distance between \eqn{cell_j} and some of its
#' nearest neighbours.
#'
#' One can choose to keep the graph as it is and consider it as a directed graph
#' (\code{do.symmetrize = FALSE}).
#' The alternative solution is to use all computed distances to extend the knn
#' graph by making the matrix symmetric. Note that connectivity graphs are
#' already symmetric, so the argument value should have no effect on the result.
#'
#' @importFrom SeuratObject DefaultAssay Cells as.Graph
#'
#' @export

setGeneric("ExpandNeighbours",
           function(object, graph.name = "RNA_nn", new.graph.name = NULL,
                    graph.type = c("distances", "connectivities"),
                    k.target = 90L, do.symmetrize = FALSE,
                    algo = c("dijkstra", "diffusion"),
                    which.dijkstra = c("auto", "igraph", "fast", "slow"),
                    dijkstra.ncores = 1L, dijkstra.tol = 1L, diffusion.iter = 26L,
                    assay = NULL, verbose = TRUE)
             standardGeneric("ExpandNeighbours"))

#' @export
#' @rdname ExpandNeighbours
setMethod("ExpandNeighbours", "Seurat",
          function(object, graph.name = "RNA_nn", new.graph.name = NULL,
                   graph.type = c("distances", "connectivities"),
                   k.target = 90L, do.symmetrize = FALSE,
                   algo = c("dijkstra", "diffusion"),
                   which.dijkstra = c("auto", "igraph", "fast", "slow"),
                   dijkstra.ncores = 1L, dijkstra.tol = 1L, diffusion.iter = 26L,
                   assay = NULL, verbose = TRUE) {
            assay <- assay %||% DefaultAssay(object)
            do.symmetrize <- do.symmetrize %||% FALSE
            graph.name_ <- paste0(graph.name, "_symmetric"[do.symmetrize])
            graph.type <- tolower(graph.type)
            graph.type <- match.arg(graph.type)
            algo <- tolower(algo)
            algo <- match.arg(algo)
            new.graph.name <- new.graph.name %||%
              sprintf("%s_%s_%dk", graph.name_, algo, k.target)
            if (algo == "dijkstra") {
              which.dijkstra <- tolower(which.dijkstra)
              which.dijkstra <- match.arg(which.dijkstra)
              expanded.mat <- expand_neighbours_dijkstra(
                object = object[[graph.name]], graph.type = graph.type,
                k.target = k.target, do.symmetrize = do.symmetrize,
                which.dijkstra = which.dijkstra, ncores = dijkstra.ncores,
                tol = dijkstra.tol, verbose = verbose
              )


            } else {
              if (! could.be.connectivity(object[[graph.name]])) {
                abort("diffusion method requires a connectivity matrix")
              }
              expanded.mat <- expand.neighbours.diffusion(
                conmat = object[[graph.name]], k.min = k.target,
                max.iter = diffusion.iter
              )
            }

            if (inherits(object[[graph.name]], "Neighbor")) {
              cells <- slot(object = object[[graph.name]], "cell.names") %||% Cells(object)
            } else {
              cells <- colnames(object[[graph.name]])
            }

            rownames(expanded.mat) <- colnames(expanded.mat) <- cells
            object[[new.graph.name]] <- as.Graph(expanded.mat)
            slot(object = object[[new.graph.name]], name = "assay.used") <- assay
            return(object)
          })

#' @keywords internal
#' @noRd
setGeneric("expand_neighbours_dijkstra",
           function(object, graph.type = c("distances", "connectivities"),
                    k.target = 90L, do.symmetrize = FALSE,
                    which.dijkstra = c("auto", "igraph", "fast", "slow"),
                    ncores = 1L, tol = 1L, verbose = TRUE)
             standardGeneric("expand_neighbours_dijkstra"))

#' @importFrom SeuratObject as.Neighbor
#' @importFrom Matrix sparseMatrix drop0
#' @keywords internal
#' @noRd
setMethod("expand_neighbours_dijkstra", "Matrix",
          function(object, graph.type = c("distances", "connectivities"),
                   k.target = 90L, do.symmetrize = FALSE,
                   which.dijkstra = c("auto", "igraph", "fast", "slow"),
                   ncores = 1L, tol = 1L, verbose = TRUE) {
            n <- ncol(object)
            const.k <- is.kconstant(object)

            vals <- unique(slot(object = object, name = "x"))
            if (length(vals) == 1 || all(vals %in% c(0,1))) {
              abort(message = paste("Provided", sQuote("Graph"),
                                    "object looks like an adjacency matrix with",
                                    "only 0s and 1s. Distances are needed"))
            }

            if (which.dijkstra == "auto") {
              if (n <= 1e4) {
                which.dijkstra <- "igraph"
              } else if (const.k) {
                which.dijkstra <- "fast"
              } else {
                which.dijkstra <- "slow"
              }
            }
            object.symmetry <- const.k.symmetry <- NULL
            igraph.mode <- "directed"
            if (do.symmetrize && !isSymmetric(object)) {
              if (which.dijkstra == "igraph") {
                igraph.mode <- "undirected"
                message('using undirected mode (igraph)\n'[verbose], appendLF = F)
              } else {
                object.symmetry <- SymmetrizeKnn(object = object)
                const.k.symmetry <- is.kconstant(object.symmetry)
                message('symmetrizing graph\n'[verbose], appendLF = F)
              }
            }

            if (graph.type == "connectivities"){
              if (! could.be.connectivity(object.symmetry %||% object)) {
                abort("provided graph doesn't seem to be a connectivity graph")
              }
              # reverse connectivities to be low when cells are near
              object <- 1 - object
              object.symmetry <- object.symmetry %iff% 1 - object.symmetry
            }

            if (which.dijkstra != "fast") {
              object <- object.symmetry %||% object
              const.k <- const.k.symmetry %||% const.k

            }

            if (which.dijkstra == "igraph") {
              if (n > 1e4) {
                warning("implementation of Dijkstra's algorithm in igraph is ",
                        "slow and not memory-efficient when too many cells/edges",
                        call. = FALSE, immediate. = TRUE)
              }
            } else if (const.k) {
              if (which.dijkstra == "slow") {
                warning("When all cells have the same number k of nearest ",
                        "neighbors, ", sQuote("fast"), " implementation of ",
                        "Dijkstra's algorithm is recommended")
              } else {
                object <- as.Neighbor(x = object)
              }
            } else if (which.dijkstra == "fast") {
              abort(message = paste(sQuote("fast"), "implementation is not",
                                    "compatible with a varying number of",
                                    "nearest neighbors per cell"))
            }
            msg <- sprintf('Computing knn (k = %d) with %s...',
                           k.target, paste("Dijkstra's algorithm using",
                                           sQuote(which.dijkstra),
                                           "implementation"))
            message(msg[verbose], appendLF = F)
            beginning <- Sys.time()
            res <- switch (which.dijkstra,
              igraph = dijkstra.igraph(knnmat = object, k.target = k.target,
                                       mode = igraph.mode, weighted = T,
                                       diag = F),
              slow = dijkstra.slow(knnmat = object, k.target = k.target,
                                   tol = tol, ncores = ncores),
              fast = dijkstra.fast(knn.neighbors = object, k.target = k.target,
                                   knn.symmetric = object.symmetry,
                                   tol = tol, ncores = ncores)
            )
            finishing <- Sys.time()
            d <- as.numeric(difftime(finishing, beginning, units = "secs"))
            dh <- floor(d / 60^2)
            dm <- floor(d / 60 - dh * 60)
            ds <- d - dh * 60^2 - dm * 60
            msg_time <- sprintf('done  |  Elapsed time: %.2d:%.2d:%s%.2f\n',
                                dh, dm, ifelse(ds < 10, '0', ''), ds)
            message(msg_time[verbose], appendLF = F)

            i <- rep(1:nrow(res$knn.idx), k.target)
            j <- as.vector(res$knn.idx)
            x <- as.vector(res$knn.dist)
            # correct igraph output when not enough neighbours
            infs <- which(is.infinite(x))
            if (length(infs) > 0) {
              x[is.infinite(x)] <- 0
              warning('Dijkstra (igraph) : could not find enough neighbours',
                      ' for ', length(infs), ' cell(s) ',
                      paste0(infs, collapse = ', '),
                      call. = F, immediate. = F)
            }
            expanded.mat <- sparseMatrix(i = i, j = j, x = x, dims = rep(n, 2))
            expanded.mat <- drop0(expanded.mat)
            if(do.symmetrize && which.dijkstra == "fast") {
              expanded.mat <- SymmetrizeKnn(expanded.mat, use.max = FALSE)
            }
            dimnames(expanded.mat) <- list(rownames(object), colnames(object))
            return(expanded.mat)
          })


#' @keywords internal
#' @noRd
setMethod("expand_neighbours_dijkstra", "Neighbor",
          function(object, graph.type = c("distances", "connectivities"),
                   k.target = 90L, do.symmetrize = FALSE,
                   which.dijkstra = c("auto", "igraph", "fast", "slow"),
                   ncores = 1L, tol = 1L, verbose = TRUE) {

            return(expand_neighbours_dijkstra(as.Graph(object),
                                              graph.type = graph.type,
                                              k.target = k.target,
                                              do.symmetrize = do.symmetrize,
                                              which.dijkstra = which.dijkstra,
                                              ncores = ncores, tol = tol,
                                              verbose = verbose))
          })

#' @importFrom Matrix rowSums diag<-
#' @keywords internal
#' @noRd
expand.neighbours.diffusion <- function(conmat, k.min, max.iter) {
  norm.factor <- 1 / max(rowSums(conmat))
  conmat <- conmat * norm.factor
  transition.mat <- cumu.transition.mat <- conmat

  k.min.found <- min(rowSums(conmat > 0))

  iter <- 2
  while(k.min.found < k.min && iter < max.iter) {
    cumu.transition.mat <- cumu.transition.mat %*% transition.mat
    conmat <- conmat + cumu.transition.mat
    iter <- iter + 1
    k.min.found <- min(rowSums(conmat > 0))
  }
  diag(conmat) <- 0
  return(conmat)
}

#' @importFrom igraph graph_from_adjacency_matrix distances
#' @keywords internal
#' @noRd
dijkstra.igraph <- function(knnmat, k.target, mode = "undirected", weighted = T, diag = F) {
  g <- graph_from_adjacency_matrix(knnmat, mode = mode, weighted = weighted,
                                   diag = diag)
  distmat <- distances(g , algorithm = "dijkstra", mode = "out")
  knn.idx  <- t(apply(distmat, 1, function(x) order(x)[1:k.target]))
  knn.dist <- t(apply(distmat, 1, function(x) x[order(x)[1:k.target]]))
  return(list(knn.idx = knn.idx, knn.dist = knn.dist))
}

# Should roots.idx & roots.dist be from the symmetrized version of the matrix ?
#' @importFrom parallel mcmapply
#' @importFrom Matrix diag<-
#' @keywords internal
#' @noRd
dijkstra.fast <- function(knn.neighbors, k.target, knn.symmetric = NULL,
                          tol = 1L, ncores = 1L) {
  knn.idx <- slot(object = knn.neighbors, name = "nn.idx")
  knn.dist <- slot(object = knn.neighbors, name = "nn.dist")
  knn.idx <- rowSort(knn.idx, knn.dist)
  knn.dist <- rowSort(knn.dist)
  fix.args <- list(knn.idx = knn.idx, knn.dist = knn.dist,
                   k=ncol(knn.idx), k.target=k.target, tol = tol)

  if (!isFALSE(knn.symmetric) && (knn.symmetric %iff% TRUE %||% FALSE)) {
    diag(knn.symmetric) <- -Inf
    i <- slot(object = knn.symmetric, name = "i") + 1
    x <- slot(object = knn.symmetric, name = "x")
    p <- slot(object = knn.symmetric, name = "p")
    j <- findInterval(seq(x)-1,p[-1]) + 1
    o <- order(i, x) # such that self is  always first
    f <- i[o]
    roots.idx  <- split(j[o], f = f)
    roots.dist <- split(x[o], f = f)
  } else {
    f <- 1:nrow(knn.idx)
    roots.idx  <- split(knn.idx,  f = f)
    roots.dist <- split(knn.dist, f = f)
  }

  if (ncores > 1) {
    res <- mcmapply(dijkstra.fast_single, roots.idx = roots.idx,
                    roots.dist = roots.dist, MoreArgs = fix.args,
                    SIMPLIFY = T, USE.NAMES = T, mc.cores = ncores)
  } else {
    res <- mapply(dijkstra.fast_single, roots.idx = roots.idx,
                  roots.dist = roots.dist, MoreArgs = fix.args,
                  SIMPLIFY = T, USE.NAMES = T)
  }
  knn.idx <- do.call(rbind, res["idx",])
  knn.dist <- do.call(rbind, res["dist",])
  return(list(knn.idx = knn.idx, knn.dist = knn.dist))
}

#' @importFrom parallel mcmapply
#' @importFrom Matrix diag<-
#' @keywords internal
#' @noRd
dijkstra.slow <- function(knnmat, k.target, tol = 1L, ncores = 1L) {
  i <- slot(object = knnmat, name = "i") + 1
  x <- slot(object = knnmat, name = "x")
  p <- slot(object = knnmat, name = "p")
  j <- findInterval(seq(x)-1, p[-1]) + 1
  v_base <- rep(Inf, ncol(knnmat))

  if (ncores > 1) {
    fix.args <- list(i = i, j = j, x = x, v = v_base,
                     k.target = k.target, tol = tol)
    res <- mcmapply(dijkstra.slow_single, 1:ncol(knnmat),
                    MoreArgs = fix.args, SIMPLIFY = TRUE, USE.NAMES = TRUE,
                    mc.cores = ncores)

  } else {
    res <- sapply(1:ncol(knnmat), dijkstra.slow_single,
                  i = i, j = j, x = x, v = v_base,
                  k.target = k.target, tol = 1L,
                  simplify = TRUE, USE.NAMES = TRUE)
  }
  knn.idx <- do.call(rbind, res["idx",])
  knn.dist <- do.call(rbind, res["dist",])
  return(list(knn.idx = knn.idx, knn.dist = knn.dist))
}

#' @keywords internal
#' @noRd
dijkstra.fast_single <- function(roots.idx, roots.dist, knn.idx, knn.dist,
                                 k, k.target, tol = 1L) {
  self <- roots.idx[1] # supposes sorted from lowest to highest distance
  reached.nodes.idx  <- roots.idx
  reached.nodes.dist <- roots.dist
  reached.nodes.dist[1] <- 0# -Inf   #current root cell always best
  roots.idx <- roots.idx[-1]
  roots.dist <- roots.dist[-1]

  # roots.idx <- roots.idx[roots.idx <= nrow(roots.idx)]
  # roots.dist <- roots.dist[!is.infinite(roots.dist)]

  no.change.since <- 0L
  tol <- max(1L, tol)
  tol.max <- tol * 50L
  iter <- 0L
  o <- order(reached.nodes.dist)[2:k.target]
  o <- o[!is.na(o)]
  reached.nodes.idx_best <- reached.nodes.idx[o]
  reached.nodes.dist_best <- reached.nodes.dist[o]
  # is.over <- (length(reached.nodes.idx_best) >= k.target && no.change.since >= tol) || no.change.since >= tol + 5
  while (length(reached.nodes.idx_best) < (k.target - 1) || no.change.since < tol) {
    iter <- iter + 1L

    if (iter > tol.max) {
      cond <- length(reached.nodes.idx_best) < (k.target - 1)
      reached.nodes.dist_best <- c(reached.nodes.dist_best,
                                   rep(0, k.target - length(reached.nodes.idx_best) - 1))
      reached.nodes.idx_best <- rep(reached.nodes.idx_best, length.out = k.target - 1)
      s1 <- c('a stable solution', 'enough neighbours')[[cond + 1]]
      s2 <- c('increase `tol` to obtain a more accurate result',
              'reduce `k.target` to obtain a constant k for all cells')[[cond + 1]]
      warn_msg <- sprintf('Dijkstra (fast) : could not find %s for cell %d (%s)',
                          s1, self, s2)
      warning(warn_msg, call. = FALSE, immediate. = FALSE)
      break
    }


    roots.nn.idx  <- knn.idx[roots.idx,]
    roots.nn.dist <- knn.dist[roots.idx,] + reached.nodes.dist[match(roots.idx, reached.nodes.idx)]
    idx <- !is.infinite(roots.nn.dist)
    roots.nn.idx  <- roots.nn.idx[idx]
    roots.nn.dist <- roots.nn.dist[idx]
    idx <- roots.nn.idx %in% reached.nodes.idx
    min.dists <- tapply(roots.nn.dist[idx], roots.nn.idx[idx], min)
    min.idx <- match(names(min.dists), reached.nodes.idx, min)
    reached.nodes.dist[min.idx] <- pmin(reached.nodes.dist[min.idx], min.dists)

    min.dists <- tapply(roots.nn.dist[!idx], roots.nn.idx[!idx], min)
    reached.nodes.dist <- c(reached.nodes.dist, unname(min.dists))
    reached.nodes.idx <- c(reached.nodes.idx, as.integer(names(min.dists)))

    o <- order(reached.nodes.dist)[2:k.target]
    o <- o[!is.na(o)]
    if (identical(reached.nodes.idx[o], reached.nodes.idx_best) &&
        identical(reached.nodes.dist[o], reached.nodes.dist_best)) {
      no.change.since <- no.change.since + 1L
    } else {
      reached.nodes.idx_best <- reached.nodes.idx[o]
      reached.nodes.dist_best <- reached.nodes.dist[o]
      no.change.since <- 0L
    }
    # idx <- which(! (reached.nodes.idx %in% roots.idx | reached.nodes.idx == self))
    roots.idx  <- reached.nodes.idx#[idx]
    roots.dist <- reached.nodes.dist#[idx]
  }
  return(list(idx = c(self, reached.nodes.idx_best),
              dist = c(0, reached.nodes.dist_best)))
}

#' @keywords internal
#' @noRd
dijkstra.slow_single <- function(self.idx, i, j, x, v, k.target, tol = 1L) {
  tol <- max(1L, tol)
  tol.max <- tol * 50L

  v[self.idx] <- 0

  idx1 <- which(i == self.idx)
  idx2 <- j[idx1]
  v[idx2] <- x[idx1]

  v_idx <- which(is.finite(v))
  v_ <- v[v_idx]
  v_o <- order(v_)
  v_o <- v_o[1:min(length(v_o), k.target)]
  reached.nodes.dist_best <- v_[v_o]
  reached.nodes.idx_best  <- v_idx[v_o]

  iter <- 0L
  no.change.since <- 0L

  while (length(reached.nodes.idx_best) < k.target || no.change.since < tol) {
    iter <- iter + 1L

    if (iter > tol.max) {
      cond <- length(reached.nodes.idx_best) < k.target
      reached.nodes.dist_best <- c(reached.nodes.dist_best,
                                   rep(0, k.target - length(reached.nodes.idx_best)))
      reached.nodes.idx_best <- rep(reached.nodes.idx_best, length.out = k.target)
      s1 <- c('a stable solution', 'enough neighbours')[[cond + 1]]
      s2 <- c('increase `tol` to obtain a more accurate result',
              'reduce `k.target` to obtain a constant k for all cells')[[cond + 1]]
      warn_msg <- sprintf('Dijkstra (slow) : could not find %s for cell %d (%s)',
                          s1, self.idx, s2)
      warning(warn_msg, call. = FALSE, immediate. = FALSE)
      break
    }

    idx3 <- which(i %in% idx2)
    idx4 <- j[idx3]
    idx4_ <- sort(unique(idx4))
    v[idx4_] <- unlist(lapply(split(pmin(v[idx4], x[idx3] + v[i[idx3]]), idx4), min),
                       recursive = F, use.names = F)

    v_idx <- which(is.finite(v))
    v_ <- v[v_idx]
    v_o <- order(v_)
    v_o <- v_o[1:min(length(v_o), k.target)]
    reached.nodes.dist_best_ <- v_[v_o]
    reached.nodes.idx_best_  <- v_idx[v_o]

    no.change.since <- no.change.since + 1L
    if (! (identical(reached.nodes.dist_best_, reached.nodes.dist_best) &&
           identical(reached.nodes.idx_best_, reached.nodes.idx_best))) {
      reached.nodes.dist_best <- reached.nodes.dist_best_
      reached.nodes.idx_best <- reached.nodes.idx_best_
      no.change.since <- 0L
    }
    idx1 <- idx3
    idx2 <- idx4_
  }
  return(list(idx = reached.nodes.idx_best,
              dist = reached.nodes.dist_best))
}
