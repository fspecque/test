#' Score a clustering result with adjusted rand index
#'
#' @description
#' Compute a score based on adjusted rand index between a clustering
#' result and one or more cell type label variable(s). 0 and 1 reflect a random
#' clustering and a perfect clustering as compared to cell type labelling
#' respectively.
#'
#' @param object A Seurat object
#' @param cell.var The name(s) of the column(s) with cell type label variable
#' (must be in the object metadata). Multiple column names are accepted
#' @param clust.var The name of the column with cluster id assignment for each
#' cell (must be in the object metadata). Only one column name is accepted
#'
#' @return \code{ScoreARI}: a named array with as many values as there are
#' common strings between cell.var and the column names of the object's
#' metadata. Names are cell.var and values are ARI.
#'
#' \code{AddScoreARI}: the updated Seurat \code{object} with the ARI score(s)
#' set for the integration.
#'
#' @export
#' @details
#' ARI is rand index corrected for chance:
#' \deqn{\displaystyle ARI = \frac{RI - RI_{expected}}{max(RI) - RI_{expected}}}
#' More precisely, a contingency table is computed with the two variables
#' \eqn{L} and \eqn{C} of \eqn{r} and \eqn{s} elements respectively. For
#' \eqn{i \in [\![1,r]\!]} and \eqn{j \in [\![1,s]\!]}, \eqn{n_{ij}} is the
#' number of common samples (i.e. cells) between \eqn{L_i} and \eqn{C_j},
#' \eqn{a_i} is the number of samples in \eqn{L_i} and \eqn{b_j} is the number
#' of samples in \eqn{C_j}. The ARI is:
#' \deqn{\displaystyle ARI = \frac{\left. \sum_{ij} \binom{n_{ij}}{2} - \left(\sum_i \binom{a_i}{2} \sum_j \binom{b_j}{2}\right) \right/ \binom{n}{2} }{ \left. \frac{1}{2} \left(\sum_i \binom{a_i}{2} + \sum_j \binom{b_j}{2}\right) - \left(\sum_i \binom{a_i}{2} \sum_j \binom{b_j}{2}\right) \right/ \binom{n}{2}}}
#'
#' @note The metric is symmetric. Switching cell.var with clust.var will return
#' the same value.
#'
#' @references Luecken, M. D., Büttner, M., Chaichoompu, K., Danese, A.,
#' Interlandi, M., Mueller, M. F., Strobl, D. C., Zappia, L., Dugas, M.,
#' Colomé-Tatché, M. & Theis, F. J. Benchmarking atlas-level data integration in
#' single-cell genomics. Nat Methods 19, 41–50 (2021).
#' \href{https://doi.org/10.1038/s41592-021-01336-8}{DOI}
#' @rdname score-ari

ScoreARI <- function(object, cell.var, clust.var = "seurat_clusters") {
  df.mtdt <- object[[]]
  if (length(clust.var) > 1) warning("Only one clustering variable allowed. Using the first one",
                                     call. = F, immediate. = T)
  clust.var <- clust.var[[1]]
  if (! clust.var %in% colnames(df.mtdt)) {
    abort(sprintf("clust.var = %s not in colnames of metadata", sQuote(clust.var)))
  }
  clust.var <- df.mtdt[, clust.var, drop = TRUE]

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

  scores <- sapply(cell.var, .ari, vec2 = clust.var,
                   simplify = "numeric", USE.NAMES = TRUE)
  return(scores)
}

#' @param integration name of the integration to score
#' @export
#' @rdname score-ari
AddScoreARI <- function(object, integration,
                        cell.var, clust.var = "seurat_clusters") {
  scores <- ScoreARI(object, cell.var = cell.var, clust.var = clust.var)

  score.names <- paste("ARI", names(scores), sep = '_')
  object <- check_misc(object)
  for (i in 1:length(scores)) {
    object <- SetMiscScore(object, integration = integration,
                           score.name = score.names[i],
                           score.value = scores[[i]],
                           class = "numeric")
  }
  return(object)
}

#' ARI: Adjusted Rand Index between two categorical variables
#' @keywords internal
#' @noRd
.ari <- function(vec1, vec2) {
  vec1 <- as.factor(vec1)
  vec2 <- as.factor(vec2)
  cross.tbl <- xtabs(~ vec1 + vec2, sparse = TRUE)

  a <- sum(choose(Matrix::rowSums(cross.tbl), 2))
  b <- sum(choose(Matrix::colSums(cross.tbl), 2))

  n <- choose(length(vec1), 2)
  ri <- sum(choose(cross.tbl@x, 2))
  ri_exp <- (a * b) / n
  ri_max <- mean(c(a,b))

  return( (ri - ri_exp) / (ri_max - ri_exp) )

}
