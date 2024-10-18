#' Score a clustering result with normalised mutual information
#'
#' @description
#' Compute a score based on normalised mutual information between a clustering
#' result and one or more cell type label variable(s). 0 and 1 reflect an absence
#' of mutual information and a perfect correlation respectively.
#'
#' @param object A Seurat object
#' @param cell.var The name(s) of the column(s) with cell type label variable
#' (must be in the object metadata). Multiple column names are accepted
#' @param clust.var The name of the column with cluster id assignment for each
#' cell (must be in the object metadata). Only one column name is accepted
#' @param average.entropy method to compute the value of the normalisation
#' denominator from each variable's entropy. one of 'mean', 'geom', 'min' and
#' 'max', namely 'arithmetic mean of', 'geometric mean of', 'minimum' and
#' 'maximum' entropy respectively.
#'
#' @return a named array with as many values as there are common strings between
#' cell.var and the column names of the object's metadata. Names are cell.var
#' and values are NMI scores.
#'
#' @export
#' @details
#' Considering a \eqn{N}-cells dataset, with \eqn{\left|L_i\right|} the number
#' of cells labelled with cell type \eqn{L_i} and \eqn{\left|C_i\right|} the
#' number of cells in cluster \eqn{C_i}. The discrete mutual information
#' \eqn{MI} approximation is given by:
#' \deqn{\displaystyle MI(L, C) = \sum_{i=1}^{\left|L\right|}\sum_{j=1}^{\left|C\right|} \left( \frac{\left|L_i \cap C_j\right|}{N} \times log \left(\frac{N \times \left|L_i \cap C_j\right|}{\left|L_i\right| \left|C_j\right|} \right) \right)}
#' Then, \eqn{MI} is normalised (scaled) by a denominator, which is computed by
#' applying a function \eqn{f} on both variables' entropies (\eqn{H}).
#' \eqn{f} can either be the arithmetic mean, geometric mean, maximum or
#' minimum of entropies.
#' \deqn{\displaystyle NMI(L, C) = \frac{MI(L, C)}{f(H(L), H(C))}}
#'
#' @note The metric is symmetric. Switching cell.var with clust.var will return
#' the same value.
#'
#' @references Luecken, M. D., Büttner, M., Chaichoompu, K., Danese, A.,
#' Interlandi, M., Mueller, M. F., Strobl, D. C., Zappia, L., Dugas, M.,
#' Colomé-Tatché, M. & Theis, F. J. Benchmarking atlas-level data integration in
#' single-cell genomics. Nat Methods 19, 41–50 (2021).
#' \href{https://doi.org/10.1038/s41592-021-01336-8}{DOI}
#'

ScoreNMI <- function(object, cell.var, clust.var = "seurat_clusters",
                     average.entropy = c('mean', 'geom', 'min', 'max')) {
  df.mtdt <- object[[]]
  if (length(clust.var) > 1) warning("Only one clustering variable allowed. Using the first one",
                                     call. = F, immediate. = T)
  clust.var <- clust.var[[1]]
  if (! clust.var %in% colnames(df.mtdt)) {
    abort(sprintf("clust.var = %s not in colnames of metadata", sQuote(clust.var)))
  }
  clust.var <- as.character(df.mtdt[, clust.var, drop = TRUE])

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
  cell.var <- asplit(df.mtdt[, cell.var, drop = FALSE], MARGIN = 2)

  scores <- sapply(cell.var, .nmi, var2 = clust.var,
                   average.entropy = average.entropy,
                   simplify = "numeric", USE.NAMES = TRUE)
  return(scores)
}

#' MI: Mutual Information between two categorical variables
#' @importFrom stats xtabs
#' @importFrom Matrix colSums rowSums
#' @keywords internal
#' @noRd
.mi <- function(vec1, vec2) {
  vec1 <- as.factor(vec1)
  vec2 <- as.factor(vec2)
  cross.tbl <- xtabs(~ vec1 + vec2, sparse = TRUE)
  i <- slot(object = cross.tbl, name = "i") + 1
  x <- slot(object = cross.tbl, name = "x")
  p <- slot(object = cross.tbl, name = "p")
  j <- findInterval(seq(x)-1, p[-1]) + 1

  n <- sum(x)
  xfreq <- x / n
  ni <- rowSums(cross.tbl)
  nj <- colSums(cross.tbl)
  outer.prod <- -log(ni[i] * nj[j]) + 2 * log(n)
  mi <- xfreq * (log(x) - log(n)) + xfreq * outer.prod
  return(max(sum(mi), 0))
}

#' Compute entropy for a categorical variable
#' @keywords internal
#' @noRd
.entropy <- function(vec) {
  counts <- tabulate(as.integer(as.factor(vec)))
  entropy <- 0
  if (length(counts) > 1) {
    n <- sum(counts)
    entropy <- -sum((counts / n) * (log(counts) - log(n)))
  }
  return(entropy)
}

#' NMI: Normalised Mutual Information between two categorical variables
#' @keywords internal
#' @noRd
.nmi <- function(var1, var2, average.entropy = c('mean', 'geom', 'min', 'max')) {
  average.entropy <- tolower(average.entropy)
  average.entropy <- match.arg(average.entropy)

  mi <- .mi(var1, var2)
  e1 <- .entropy(var1)
  e2 <- .entropy(var2)

  denorm <- switch(
    average.entropy,
    mean = mean(c(e1, e2)),
    geom = sqrt(e1 * e2),
    min = min(e1, e2),
    max = max(e1, e2))

  return( mi / denorm )
}
