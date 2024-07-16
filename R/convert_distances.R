#' @include utils.R
NULL

#' Smooth knn distances
#'
#' @description
#' Directly adapted from \pkg{UMAP}'s function \code{smooth_knn_dist()}.
#'
#' \emph{Compute a continuous version of the distance to the kth nearest neighbor.
#' That is, this is similar to knn-distance but allows continuous k values rather
#' than requiring an integral k. In essence we are simply computing the distance
#' such that the cardinality of fuzzy set we generate is k.}
#'
#' @param distances matrix of distances of k columns and a row per cell
#' @param k desired number of nearest neighbours to approximate.
#' If \code{NULL}, set to number of columns of matrix \code{distances} (default)
#' @param max.iter maximum number of iterations during smoothing
#' @param local.connectivity The local connectivity required, i.e. the expected
#' number of nearest neighbours locally connected. The higher, the more local
#' connections are output. In practice, this should be not more than the local
#' intrinsic dimension of the manifold.
#' @param bandwidth The target bandwidth of the kernel. Larger values will
#' produce larger return values.
#'
#' @return A list containing 2 vectors with one value per cell:
#' \itemize{
#'   \item \strong{sigmas}: approximated distance to 1st nearest neighbor of
#'   each cell.
#'   \item \strong{rhos}: distance to the kth nearest neighbor of each cell.
#' }
#'
#' @note This function is adapted from \pkg{umap-learn} original python
#' implementation. You can have a look at the
#' \href{https://umap-learn.readthedocs.io/en/latest/api.html#umap.umap_.smooth_knn_dist}{original documentation}
#'
#' @references McInnes, L., Healy, J. & Melville, J. UMAP: Uniform Manifold
#' Approximation and Projection for Dimension Reduction. arXiv preprint (2018).
#' \href{https://doi.org/10.48550/ARXIV.1802.03426}{DOI}
#'
#' @seealso \code{\link{compute_membership_strengths}} and
#' \code{\link{fuzzy_simplicial_set}}
#' @keywords internal
#' @noRd

smooth_knn_dist <- function(distances, k = NULL, max.iter=64,
                            local.connectivity=1.0, bandwidth=1.0) {
  n <- nrow(distances)
  k <- k %||% ncol(distances)
  N <- n * k

  SMOOTH_K_TOLERANCE <- 1e-5
  MIN_K_DIST_SCALE <- 1e-3

  target = log2(k) * bandwidth
  rho <- rep(0, n)
  mean.distance <- mean(distances)
  mean.distances <- rowMeans(distances)

  non.zero.dists <- lapply(split(distances, 1:n), function(v) v[v > 0])
  n.non.zeros <- rowSums(distances > 0)
  idx <- n.non.zeros >= local.connectivity
  index <- as.integer(floor(local.connectivity))
  interpolation <- local.connectivity - index
  cond1 <- index > 0
  cond2 <- interpolation > SMOOTH_K_TOLERANCE

  if (any(idx)) {
    if (cond1) {
      v1 <- sapply(non.zero.dists[idx], function(v) v[index],
                   simplify = "array", USE.NAMES = FALSE)
      rho[idx] <- v1
      if (cond2) {
        v2 <- sapply(non.zero.dists[idx], function(v) v[index + 1],
                     simplify = "array", USE.NAMES = FALSE)
        rho[idx] <- rho[idx] + interpolation * (v2 - v1)
      }
    } else {
      v3 <- sapply(non.zero.dists[idx], function(v) v[1],
                   simplify = "array", USE.NAMES = FALSE)
      rho[idx] <- interpolation * v3
    }
  }
  idx <- !idx & n.non.zeros > 0
  if (any(idx)) {
    rho[idx] <- sapply(non.zero.dists[idx], max,
                       simplify = "array", USE.NAMES = FALSE)
  }

  idx_continue <- rep(T, n)
  lo <- rep(0.0, n)
  hi <- lo + Inf
  mid <- lo + 1
  distances <- distances - rho
  for (it in 1:max.iter) {
    psum <- rep(0.0, n)
    for(j in 2:k) {
      d <- distances[, j]
      idx_ <- d > 0
      idx <- idx_ & idx_continue
      if (any(idx)) {
        psum[idx] <- psum[idx] + exp(-(d[idx] / mid[idx]))
      }
      idx <- !idx_ & idx_continue
      if (any(idx)) {
        psum[idx] <- psum[idx] + 1
      }
    }

    idx_continue[idx_continue] <- abs(psum[idx_continue] - target) >= SMOOTH_K_TOLERANCE

    if (! any(idx_continue)) {
      break
    }

    idx_ <- psum > target
    idx <- idx_ & idx_continue
    hi[idx] <- mid[idx]
    mid[idx] <- (lo[idx] + hi[idx]) / 2.0
    idx <- !idx_ & idx_continue
    lo[idx] <- mid[idx]
    iidx <- idx & is.infinite(hi)
    mid[iidx] <- mid[iidx] * 2
    iidx <- idx & !is.infinite(hi)
    mid[iidx] <- (lo[iidx] + hi[iidx]) / 2.0
  }
  result <- mid
  idx_ <- rho > 0
  scaled.mean.dist  <- MIN_K_DIST_SCALE * mean.distance
  scaled.mean.dists <- MIN_K_DIST_SCALE * mean.distances
  idx <- idx_ & result < scaled.mean.dists
  result[idx] <- scaled.mean.dists[idx]
  idx <- !idx_ & result < scaled.mean.dist
  result[idx] <- scaled.mean.dist

  return(list(sigmas = result, rhos = rho))
}


#' Convert knn distances to connectivities
#'
#' @description
#' Directly adapted from \pkg{UMAP}'s function \code{compute_membership_strengths()}.
#'
#' \emph{Construct the membership strength data for the 1-skeleton of each local
#' fuzzy simplicial set – this is formed as a sparse matrix where each row is a
#' local fuzzy simplicial set, with a membership strength for the 1-simplex to
#' each other data point.}
#'
#' @param knn.idx matrix of indices of nearest neighbours of k columns and a row
#' per cell, corresponding to \code{knn.dist}
#' @param knn.dist matrix of distances of k columns and a row per cell
#' @param sigmas normalization factor derived from the metric tensor
#' approximation. Corresponds to \code{smooth_knn_dist()$sigmas}
#' @param rhos local connectivity adjustment. Corresponds to
#' \code{smooth_knn_dist()$rhos}
#' @param return.dists whether to return the matrix of distances (\code{FALSE}
#' by default)
#' @param bipartite whether the knn network is bipartite (\code{FALSE} by
#' default)
#'
#' @return A list containing 3 vectors (or 4 when \code{return.dists = TRUE})
#' with one value per cell:
#' \itemize{
#'   \item \strong{rows}: row indices of the resulting sparse matrix (i)
#'   \item \strong{cols}: column indices of the resulting sparse matrix (j)
#'   \item \strong{vals}: values of the resulting sparse matrix for every (i,j)
#'   coordinates
#'   \item \strong{dists}: \code{knn.dist} as a vector
#' }
#'
#' @note This function is adapted from \pkg{umap-learn} original python
#' implementation. You can have a look at the
#' \href{https://umap-learn.readthedocs.io/en/latest/api.html#umap.umap_.compute_membership_strengths}{original documentation}
#'
#' @references McInnes, L., Healy, J. & Melville, J. UMAP: Uniform Manifold
#' Approximation and Projection for Dimension Reduction. arXiv preprint (2018).
#' \href{https://doi.org/10.48550/ARXIV.1802.03426}{DOI}
#'
#' @seealso \code{\link{smooth_knn_dist}} and \code{\link{fuzzy_simplicial_set}}
#' @keywords internal
#' @noRd

compute_membership_strengths <- function(knn.idx, knn.dist, sigmas, rhos,
                                         return.dists=FALSE, bipartite=FALSE) {
  n <- nrow(knn.idx)
  k <- ncol(knn.idx)
  N <- n * k

  rows <- rep(1:n, k)
  cols <- as.vector(knn.idx)
  dists <-  NULL

  vals <- rep(0, N)
  idx.not.self <- rep(T, N)
  if (!bipartite) {
    idx.not.self <- knn.idx != 1:n
  }
  idx <- knn.dist - rhos <= 0 | sigmas == 0
  vals[idx.not.self &  idx] <- 1
  vals[idx.not.self & !idx] <- exp(-((knn.dist - rhos) / sigmas))[idx.not.self & !idx]
  if (return.dists) {
    dists <- as.vector(knn.dist)
  }
  return(list(rows = rows, cols = cols, vals = vals, distances = dists))
}


#' Compute a fuzzy simplicial set and performs a fuzzy union
#'
#' @description
#' Directly adapted from \pkg{UMAP}'s function \code{fuzzy_simplicial_set()}.
#'
#' \emph{Given [...] a neighborhood size, and a measure of distance compute the
#' fuzzy simplicial set (here represented as a fuzzy graph in the form of a
#' sparse matrix) associated to the data. This is done by locally approximating
#' geodesic distance at each point, creating a fuzzy simplicial set for each
#' such point, and then combining all the local fuzzy simplicial sets into a
#' global one via a fuzzy union.}
#'
#' @param knn.idx matrix of indices of nearest neighbours of k columns and a row
#' per cell, corresponding to \code{knn.dist}. \strong{Must be sorted in a
#' row-wise manner by \code{knn.dist}}
#' @param knn.dist matrix of distances of k columns and a row per cell.
#' \strong{Must be sorted in a row-wise manner by \code{knn.dist}}
#' @param set.op.mix.ratio Interpolate between (fuzzy) union and intersection as
#' the set operation used to combine local fuzzy simplicial sets to obtain a
#' global fuzzy simplicial set. Both fuzzy set operations use the product
#' t-norm. The value of this parameter should be between 0.0 and 1.0; a value of
#' 1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy intersection.
#' @param local.connectivity The local connectivity required, i.e. the expected
#' number of nearest neighbours locally connected. The higher, the more local
#' connections are output. In practice, this should be not more than the local
#' intrinsic dimension of the manifold.
#' @param niter.smoothing maximum number of iterations during smoothing
#' @param apply.set.operations whether to perform a fuzzy union or intersection
#' of local fuzzy simplicial sets into a global fuzzy simplicial set.
#' \code{TRUE} by default
#' @param bipartite whether the knn network is bipartite (\code{FALSE} by
#' default)
#' @param return.dists whether to return the matrix of distances (\code{FALSE}
#' by default)
#' @param verbose whether to print progress messages
#'
#' @return A list containing 3 vectors (or 4 when \code{return.dists = TRUE})
#' with one value per cell:
#' \itemize{
#'   \item \strong{connectivities}: sparse matrix of computed connectivities of
#'   size ncells x ncells.
#'   \item \strong{sigmas}: approximated distance to 1st nearest neighbor of
#'   each cell.
#'   \item \strong{rhos}: distance to the kth nearest neighbor of each cell.
#'   \item \strong{distances}: sparse matrix of distances (from \code{knn.dist})
#'   of size ncells x ncells.
#' }
#'
#' @importFrom Matrix t sparseMatrix drop0
#'
#' @note This function is adapted from \pkg{umap-learn} original python
#' implementation. You can have a look at the
#' \href{https://umap-learn.readthedocs.io/en/latest/api.html#umap.umap_.fuzzy_simplicial_set}{original documentation}
#'
#' @references McInnes, L., Healy, J. & Melville, J. UMAP: Uniform Manifold
#' Approximation and Projection for Dimension Reduction. arXiv preprint (2018).
#' \href{https://doi.org/10.48550/ARXIV.1802.03426}{DOI}
#' @keywords internal
#' @noRd

# /!\ knn.dist & knn.idx must be sorted by knn.dist
fuzzy_simplicial_set <- function(knn.dist, knn.idx, set.op.mix.ratio = 1,
                                 local.connectivity = 1, niter.smoothing = 64L,
                                 apply.set.operations = TRUE, bipartite = FALSE,
                                 return.dists = NULL, verbose = TRUE) {
  k <- ncol(knn.dist)
  n <- nrow(knn.dist)
  return.dists <- (return.dists %iff% as.logical(return.dists)) %||% FALSE

  e <- environment()

  message("1) Smoothing distances\n"[verbose], appendLF = F)
  list2env(smooth_knn_dist(distances = knn.dist, k = k, max.iter = niter.smoothing,
                           local.connectivity = local.connectivity), envir = e)

  message("2) Computing connectivities\n"[verbose], appendLF = F)
  list2env(compute_membership_strengths(knn.idx = knn.idx, knn.dist = knn.dist,
                                        sigmas = sigmas, rhos = rhos,
                                        return.dists = return.dists,
                                        bipartite = bipartite), envir = e)

  res.mat <- drop0(sparseMatrix(i=rows, j=cols, x=vals, dims = rep(n, 2)))

  if (apply.set.operations) {
    res.maT <- t(res.mat)
    prod.mat <- res.mat * res.maT
    res.mat <- (set.op.mix.ratio * (res.mat + res.maT - prod.mat) +
                  (1.0 - set.op.mix.ratio) * prod.mat)
  }

  res.l <- list(connectivities = res.mat, sigmas = sigmas, rhos = rhos,
                distances = NULL)
                # distances = list(rows = rows, cols = cols, x=distances))

  if (return.dists) {
    message("3) Computing distances\n"[verbose], appendLF = F)
    res.l$distances <- symmetrize.pmax.sparse(i=rows, j=cols, x=distances,
                                              height = n)
  }

  return(res.l)
}


#' Wrapper of \code{fuzzy_simplicial_set()}
#'
#' @description
#' Wraps \code{fuzzy_simplicial_set()} to only output the sparse matrix of
#' connectivities and to ensure that the distances are adequately sorted and
#' that matrices of distances and indices correspond. This is an alternative to
#' the Gaussian kernel method from \pkg{scanpy}.
#'
#' @param knn.idx matrix of indices of nearest neighbours of k columns and a row
#' per cell, corresponding to \code{knn.dist}.
#' @param knn.dist matrix of distances of k columns and a row per cell.
#' @param sorted.dist whether \code{knn.dist} and \code{knn.idx} are sorted in a
#' row-wise manner by \code{knn.dist}. Matrices are sorted by default, should
#' only be disabled when matrices have already been sorted.
#' @param set.op.mix.ratio Interpolate between (fuzzy) union and intersection as
#' the set operation used to combine local fuzzy simplicial sets to obtain a
#' global fuzzy simplicial set. Both fuzzy set operations use the product
#' t-norm. The value of this parameter should be between 0.0 and 1.0; a value of
#' 1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy intersection.
#' @param local.connectivity The local connectivity required, i.e. the expected
#' number of nearest neighbours locally connected. The higher, the more local
#' connections are output. In practice, this should be not more than the local
#' intrinsic dimension of the manifold.
#' @param niter.smoothing maximum number of iterations during smoothing
#' @param apply.set.operations whether to perform a fuzzy union or intersection
#' of the local fuzzy simplicial sets into a global fuzzy simplicial set.
#' \code{TRUE} by default
#' @param bipartite whether the knn network is bipartite (\code{FALSE} by
#' default)
#' @param return.dists whether to return the matrix of distances (\code{FALSE}
#' by default)
#' @param verbose whether to print progress messages
#'
#' @return A sparse matrix of computed connectivities of size ncells x ncells.
#'
#' @seealso \code{\link{compute.gauss.connectivities}}
#' @keywords internal
#' @noRd

compute.umap.connectivities <- function(knn.dist, knn.idx, sorted.dist = FALSE,
                                        set.op.mix.ratio = 1,
                                        local.connectivity = 1,
                                        niter.smoothing = 64L,
                                        apply.set.operations = TRUE,
                                        bipartite = FALSE,
                                        return.dists = NULL, verbose = TRUE) {
  sort.by.dist <- !isTRUE(sorted.dist %||% FALSE)
  k <- unique(ncol(knn.dist), ncol(knn.idx))
  n <- unique(nrow(knn.dist), nrow(knn.idx))

  msg.discrepancy <- sprintf("%s does not correspond to %s",
                             sQuote("knn.dist"), sQuote("knn.idx"))
  if (length(k) > 1) {
    abort(message = paste(msg.discrepancy, "(different k values)"))
  }
  if (length(n) > 1) {
    abort(message = paste(msg.discrepancy, "(different number of cells)"))
  }

  if (sort.by.dist) {
    knn.idx  <- rowSort(knn.idx, by = knn.dist, ncol = k)
    knn.dist <- rowSort(knn.idx, ncol = k)
  }
  res <- fuzzy_simplicial_set(knn.dist, knn.idx, set.op.mix.ratio = set.op.mix.ratio,
                              local.connectivity = local.connectivity,
                              niter.smoothing = niter.smoothing,
                              apply.set.operations = apply.set.operations,
                              bipartite = bipartite, return.dists = return.dists,
                              verbose = verbose)
  return(res$connectivities)
}


#' Compute a Gaussian kernel to estimate connectivities from distances
#'
#' @description
#' Directly adapted from \pkg{scanpy}'s function \code{gauss()} to compute
#' connectivities. This is an alternative to fuzzy simplicial set method from
#' \pkg{UMAP}.
#'
#' \emph{Derive Gaussian connectivities between data points from their distances.}
#'
#' @param knn.idx matrix of indices of nearest neighbours of k columns and a row
#' per cell, corresponding to \code{knn.dist}.
#' @param knn.dist matrix of distances of k columns and a row per cell.
#' @param sorted.dist whether \code{knn.dist} and \code{knn.idx} are sorted in a
#' row-wise manner by \code{knn.dist}. Matrices are sorted by default, should
#' only be disabled when matrices have already been sorted.
#' @param sigmas by default, a sigma value (i.e. width of the kernel) per cell
#' is computed internally. Those sigmas are controlling for each cell's
#' connectivities range and magnitude. Alternatively, you can provide your own
#' width(s). If the length of \code{sigmas} is shorter than the number of cells,
#' it will be recycled.
#' @param median.sigma which estimation method to use when \code{sigmas = NULL}.
#' By default, use the cell’s distance to its kth nearest neighbor. When
#'  \code{median.sigma = TRUE}, use the median of distances to the nearest
#'  neighbors (excluding self)
#' @param verbose whether to print progress messages
#'
#' @return A sparse matrix of computed connectivities of size ncells x ncells.
#'
#' @note This function is adapted from \pkg{scanpy} original python
#' implementation. You can have a look at the
#' \href{https://github.com/scverse/scanpy/blob/b918a23eb77462837df90d7b3a30a573989d4d48/src/scanpy/neighbors/_connectivity.py#L18}{original function}\cr
#' \code{median.sigmas = FALSE} is the method used in Haghverdi L. et al.,2016.
#'
#' @references Hie, B., Bryson, B. & Berger, B. Efficient integration of
#' heterogeneous single-cell transcriptomes using Scanorama. Nat Biotechnol 37,
#' 685–691 (2019). \href{https://doi.org/10.1038/s41587-019-0113-3}{DOI}
#'
#' @seealso \code{\link{compute.umap.connectivities}}
#' @keywords internal
#' @noRd

compute.gauss.connectivities <- function(knn.dist, knn.idx, sorted.dist = FALSE,
                                         sigmas=NULL, median.sigma = FALSE) {
  sort.by.dist <- !isTRUE(sorted.dist %||% FALSE)
  k <- unique(ncol(knn.dist), ncol(knn.idx))
  n <- unique(nrow(knn.dist), nrow(knn.idx))

  msg.discrepancy <- sprintf("%s does not correspond to %s",
                             sQuote("knn.dist"), sQuote("knn.idx"))
  if (length(k) > 1) {
    abort(message = paste(msg.discrepancy, "(different k values)"))
  }
  if (length(n) > 1) {
    abort(message = paste(msg.discrepancy, "(different number of cells)"))
  }

  if (sort.by.dist) {
    knn.idx  <- rowSort(knn.idx, by = knn.dist, ncol = k)
    knn.dist <- rowSort(knn.idx, ncol = k)
  }

  sigmas <- sigmas %iff% if (isFALSE(sigmas)){ NULL } %||% {
    if (median.sigma) {
      # knn.dist is sorted, we remove the nearest neighbour per cell (i.e. self)
      apply(knn.dist[,-1], 1, median)
    } else {
      # knn.dist is sorted, we take the distance of the furthest neighbour per cell
      knn.dist[, k, drop = TRUE] / 2
    }
  }

  if (length(sigmas) != n) {
    warn.msg <- c(paste("more sigmas than cells, taking the first", n),
                  "number of cells is not a multiple of sigmas length, recycling sigmas")
    warning(warn.msg[(length(sigmas) < n) + 1])
    sigmas <- rep(sigmas, length.out = n)
  }
  sigmas2 <- sigmas^2

  bw.approx <- sigmas2 + sigmas2[knn.idx]
  bw.exact  <- sigmas  * sigmas[knn.idx] * 2

  correct.factor <- sqrt(bw.exact / bw.approx)
  conns <- correct.factor * exp(-knn.dist^2 / bw.approx)

  return(
    symmetrize.pmax.sparse(i = rep(1:n, times = k),
                           j = as.vector(knn.idx),
                           x = as.vector(conns),
                           height = nrow(conns))
  )
}

#' Derive connectivities from distances
#'
#' @description
#' Adapted from \pkg{scanpy}'s strategy in \code{sc.pp.neighbors()} to compute
#' connectivities. Two methods are available, using a Gaussian kernel or a fuzzy
#' union of simplical sets as implemented in \pkg{umap-learn}
#'
#' @param object a \code{Seurat} object
#' @param neighbors.use name of a \code{Neighbor} instance stored in the
#' \code{Seurat object} to derive connectivities from.
#' @param graph.name name of the return \code{Graph} of connectivities to store
#' in the \code{Seurat object}.
#' @param assay name of the assay to reference in the output \code{Graph} object.
#' Use the default assay of \code{object} if not provided.
#' @param umap.set.op.mix.ratio float between 0 and 1. Controls how the fuzzy
#' sets are mixed to obtain a global fuzzy simplicial set. 0 and 1 correspond to
#' a pure fuzzy intersection and union, respectively. Both fuzzy set operations
#' use the product t-norm. Only applies to the \strong{umap} method.
#' @param umap.local.connectivity The local connectivity required, i.e. the
#' expected number of nearest neighbours locally connected. The higher, the more
#' local connections are output. In practice, this should be not more than the
#' local intrinsic dimension of the manifold. Only applies to the \strong{umap}
#' method when \code{umap.apply.set.operations = TRUE}.
#' @param umap.niter.smoothing maximum number of iterations during the smoothing
#' process of distances. Only applies to the \strong{umap} method.
#' @param umap.apply.set.operations set to \code{FALSE} to disable the fuzzy
#' union or intersection of the local fuzzy simplicial sets into a global fuzzy
#' simplicial set. Only applies to the \strong{umap}
#' method.
#' @param umap.bipartite whether the knn network is bipartite (\code{FALSE} by
#' default). Only applies to the \strong{umap} method.
#' @param gauss.sigmas by default, a sigma value (i.e. width of the kernel) per
#' cell is computed internally. Those sigmas are controlling for each cell's
#' connectivities range and magnitude. Alternatively, you can provide your own
#' width(s). If the length of \code{sigmas} is shorter than the number of cells,
#' it will be recycled. Only applies to the \strong{Gaussian} method.
#' @param gauss.median.sigma which estimation method to use when
#' \code{sigmas = NULL}. By default, use the cell’s distance to its kth nearest
#' neighbour. When \code{median.sigma = TRUE}, use the median of distances to
#' the nearest neighbours (excluding self). Only applies to the
#' \strong{Gaussian} method.
#' @param verbose whether to print progress messages
#'
#' @return the Seurat object with a new \code{Graph} instance of name
#' \code{graph.name}
#'
#' @export
#'
#' @importFrom SeuratObject DefaultAssay Cells as.Graph
#' @note The UMAP method is a re-implementation of the function
#' \code{fuzzy_simplicial_set()} from \pkg{umap-learn} and should estimate
#' identical connectivities. You can check the
#' \href{https://umap-learn.readthedocs.io/en/latest/api.html#umap.umap_.fuzzy_simplicial_set}{original documentation}.\cr\cr
#' The Gaussian kernel method is a re-implementation of the analogous function
#' from \pkg{scanpy} called by \code{\code{sc.pp.neighbors()}}. You can have a
#' look at the
#' \href{https://github.com/scverse/scanpy/blob/b918a23eb77462837df90d7b3a30a573989d4d48/src/scanpy/neighbors/_connectivity.py#L18}{original function}\cr
#' \code{median.sigmas = FALSE} is the method used in Haghverdi L. et al.,2016.
#'
#' @references McInnes, L., Healy, J. & Melville, J. UMAP: Uniform Manifold
#' Approximation and Projection for Dimension Reduction. arXiv preprint (2018).
#' \href{https://doi.org/10.48550/ARXIV.1802.03426}{DOI}
#' @references Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale
#' single-cell gene expression data analysis. Genome Biology 19, (2018).
#' @references Coifman, R. R., Lafon, S., Lee, A. B., Maggioni, M., Nadler, B.,
#' Warner, F. & Zucker, S. W. Geometric diffusions as a tool for harmonic
#' analysis and structure definition of data: Diffusion maps. PNAS 102, 7426–7431
#' (2005). \href{https://doi.org/10.1073/pnas.0500334102}{DOI}
#' @references Haghverdi, L., Büttner, M., Wolf, F. A., Buettner, F. & Theis, F.
#' J. Diffusion pseudotime robustly reconstructs lineage branching. Nature
#' Methods 13, 845–848 (2016). \href{https://doi.org/10.1038/nmeth.3971}{DOI}

GetConnectivities <- function(object, neighbors.use, method = c("umap", "gauss"),
                              graph.name = NULL, assay = NULL,
                              umap.set.op.mix.ratio = 1,
                              umap.local.connectivity = 1,
                              umap.niter.smoothing = 64L,
                              umap.apply.set.operations = TRUE,
                              umap.bipartite = FALSE,
                              gauss.sigmas=NULL, gauss.median.sigma = FALSE,
                              verbose = TRUE) {
  method <- tolower(method)
  method <- match.arg(method)
  method.msg <- switch (method,
    "umap" = "UMAP's fuzzy simplicial set",
    "gauss" = "a gaussian kernel"
  )
  graph.name <- graph.name %||% paste0("connectivities_", method, "_", neighbors.use)
  assay <- assay %||% DefaultAssay(so)

  graph <- object[[neighbors.use]]
  if (class(graph) != "Neighbor") {
    rlang::abort(sprintf("%s object required, got a %s object instead",
                         sQuote("Neighbor"), sQuote(class(graph))))
  }

  knn.idx  <- slot(object = graph, name = "nn.idx")
  knn.dist <- slot(object = graph, name = "nn.dist")
  dist.metric <- slot(object = graph, name = "alg.info")$metric
  sorted.dist <- all(rowSorted(knn.dist))

  message(paste("Computing connectivities from", dist.metric, "distances using",
          method.msg, "\n")[verbose], appendLF = FALSE)

  conns <- switch (method,
    "umap" = compute.umap.connectivities(knn.dist = knn.dist, knn.idx = knn.idx,
                                         sorted.dist = sorted.dist,
                                         set.op.mix.ratio = umap.set.op.mix.ratio,
                                         local.connectivity = umap.local.connectivity,
                                         niter.smoothing = umap.niter.smoothing,
                                         apply.set.operations = umap.apply.set.operations,
                                         bipartite = umap.bipartite,
                                         return.dists = FALSE, verbose = verbose),
    "gauss" = compute.gauss.connectivities(knn.dist = knn.dist, knn.idx = knn.idx,
                                           sorted.dist = sorted.dist,
                                           sigmas = gauss.sigmas,
                                           median.sigma = gauss.median.sigma)
  )
  colnames(conns) <- rownames(conns) <- Cells(so)
  object[[graph.name]] <- as.Graph(conns)
  slot(object = object[[graph.name]], name = "assay.used") <- assay
  message(paste("Connectivity Graph saved as", sQuote(graph.name), "\n")[verbose],
          appendLF = FALSE)
  return(object)
}
