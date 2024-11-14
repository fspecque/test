#' @include utils.R
NULL
# \deqn{score = \sum_{c \in L} \left( \frac{1}{\left| E \right| + \left| E\prime \right|} \times \frac{\left| E \right|}{\left| E\prime \right|}\right)}

#' Score a knn graph based on cell-type label connectivity
#'
#' @description
#' Compute scores based on the connectivity between cells sharing the same label.
#' The score can be calculated in two fashions, either by quantifying the
#' connectivity of the largest subgraph for each cell label (identical to the
#' score used in Luecken M.D. \emph{et al.}, 2022), or directly on the whole
#' graph.
#'
#' @param object A Seurat object
#' @param graph.name The name of the knn graph to score.
#' @param cell.var The name of the  cell variable (must be in the object metadata).
#' @param do.symmetrize whether to symmetrize the knn graph. Set to\code{FALSE}
#' to disable (not recommended, see \strong{Details})
#' @param per.component whether to use the same score as in Luecken M.D.
#' \emph{et al.}, 2022. \code{TRUE} by default. See \strong{Details}.
#' @param count.self whether to account for loops (i.e. cells connected to
#' themselves). \code{FALSE} by default
#' @param weight.by.ncells whether to weight the connectivity-derived scores
#' computed on each cell type label by their relative proportion in number of
#' cells. By default (\code{FALSE}), the overall score is computed as the mean
#' of scores computed per label. Ignored when \code{per.component = TRUE}
#'
#' @return \code{ScoreConnectivity}: a single float between 0 and 1,
#' corresponding to the connectivity score.
#'
#' \code{AddScoreConnectivity}: the updated Seurat \code{object} with the Graph
#' connectivity score set for the integration.
#'
#' @importFrom SeuratObject Neighbors Graphs
#'
#' @export
#' @details
#' The default parameters (\code{per.component = TRUE, count.self = FALSE})
#' correspond to the original score from Luecken M.D. \emph{et al.}, 2022. It is
#' computed as follow: for each cell type label \eqn{c} among the set of all
#' labels \eqn{L}, the sub graph \eqn{subG_c} exclusively composed of cells
#' \eqn{c} is extracted from the full graph \eqn{G}. Then, the size (i.e. the
#' number of cells, hence of vertices \eqn{V}) of the largest connected
#' component is divided by the size of \eqn{subG_c}. Then, the mean of sub
#' graphs' scores is computed:
#' \deqn{ratio_c = \frac{max(\left|V(CC(subG_c))\right|)}{\left|V(subG_c)\right|}}
#' \deqn{score = \frac{1}{\left| L \right|}\sum_{c \in L} ratio_c}
#'
#'
#' When \code{per.component = FALSE}, the connectivity is computed on the full
#' graph \eqn{G}. Let's consider the set of all labels \eqn{L},
#' \eqn{c \in L} and \eqn{L\prime = L \setminus \{c\}}. Let's also denote the
#' edges between cells \eqn{\in \{c\}} \eqn{E_c = E_{c \to c}(G)}
#' and the edges connecting cells \eqn{\in \{c\}} with cells \eqn{\in L\prime}
#' \eqn{E_c^{\prime} = E_{c \to L\prime}(G)}.
#' When \code{weight.by.ncells = TRUE}, the score is computed as follow:
#' \deqn{score = \sum_{c \in L} \left( \frac{\left| V(subG_c) \right|}{\left| V(G) \right|} \times \frac{\left| E_c \right|}{\left| E_c^\prime \right| + \left| E_c \right|}\right)}
#'
#' When \code{weight.by.ncells = FALSE}, the score is the mean of ratio of edges:
#' \deqn{score = \frac{1}{\left| L \right|}\sum_{c \in L} \frac{\left| E_c \right|}{\left| E_c^\prime \right| + \left| E_c \right|}}
#'
#' In either case, it is recommended to keep \code{do.symmetrize = TRUE}.
#'
#' @note This score is an adaptation of the graph connectivity score as
#' described in Luecken M.D. \emph{et al.}, 2022.
#'
#' @references Luecken, M. D., Büttner, M., Chaichoompu, K., Danese, A.,
#' Interlandi, M., Mueller, M. F., Strobl, D. C., Zappia, L., Dugas, M.,
#' Colomé-Tatché, M. & Theis, F. J. Benchmarking atlas-level data integration in
#' single-cell genomics. Nat Methods 19, 41–50 (2021).
#' \href{https://doi.org/10.1038/s41592-021-01336-8}{DOI}
#'
#' @seealso \code{\link[Seurat]{FindNeighbors}}
#' @rdname score-connectivity

ScoreConnectivity <- function(object, graph.name, cell.var, do.symmetrize = TRUE,
                              per.component = TRUE, count.self = FALSE,
                              weight.by.ncells = FALSE) {
  if (! graph.name %in% c(Graphs(object), Neighbors(object))) {
    rlang::abort(paste(sQuote(graph.name), "not in the Seurat object"))
  }
  if (!cell.var %in% colnames(object[[]])) {
    abort(paste(sQuote(cell.var), "not in the Seurat object's metadata"))
  }
  cell.var <- object[[]][, cell.var, drop = FALSE]

  scores <- .graph.connectivity(object[[graph.name]], cell.var = cell.var,
                                do.symmetrize = do.symmetrize,
                                per.component = per.component,
                                count.self = count.self)
  if (weight.by.ncells && ! per.component) {
    cells.prop <- proportions(table(cell.var[, 1, drop = TRUE]))[names(scores)]
    scores <- sum(scores * cell.props)
  } else {
    scores <- mean(scores)
  }
  return(scores)
}

#' @param integration name of the integration to score
#' @export
#' @rdname score-connectivity
AddScoreConnectivity <- function(object, integration,
                                 graph.name, cell.var, do.symmetrize = TRUE,
                                 per.component = TRUE, count.self = FALSE,
                                 weight.by.ncells = FALSE) {
  scores <- ScoreConnectivity(object, graph.name = graph.name,
                              cell.var = cell.var,
                              do.symmetrize = do.symmetrize,
                              per.component = per.component,
                              count.self = count.self,
                              weight.by.ncells = weight.by.ncells)

  object <- check_misc(object)
  object <- SetMiscScore(object, integration = integration,
                         score.name = "Graph.connectivty",
                         score.value = scores)
  return(object)
}

#' @keywords internal
setGeneric(".graph.connectivity",
           function(object, cell.var, do.symmetrize, per.component, count.self)
           standardGeneric(".graph.connectivity"))


#' @importFrom SeuratObject as.Graph
#' @keywords internal
setMethod(".graph.connectivity", "Neighbor",
          function(object, cell.var, do.symmetrize, per.component, count.self) {
            g <- as.Graph(object)
            .graph.connectivity(g, cell.var = cell.var,
                                do.symmetrize = do.symmetrize,
                                per.component = per.component,
                                count.self = count.self)
          })

#' @importFrom Matrix isSymmetric diag rowSums
#' @importFrom igraph components graph_from_adjacency_matrix
#' @keywords internal
setMethod(".graph.connectivity", "Matrix",
          function(object, cell.var, do.symmetrize, per.component, count.self) {
            object <- object[rownames(cell.var), rownames(cell.var)]
            diag(object) <- as.integer(count.self)
            if (do.symmetrize && !isSymmetric(object)) {
              object <- SymmetrizeKnn(object, use.max = FALSE)
            }
            object <- object > 0
            if (per.component) {
              cell.var <- as.character(cell.var[, 1, drop = TRUE])
              scores <- sapply(unique(cell.var), function(cell) {
                cells <- which(cell.var == cell)
                csize <- components(
                  graph_from_adjacency_matrix(object[cells, cells],
                                              weighted = NULL),
                  mode = 'weak')$csize
                max(csize) / sum(csize)
              })
            } else {
              cnts <- GetNeighborsPerBatch(object, cell.var, count.self = count.self)
              scores <- diag(cnts) / rowSums(cnts)
              names(scores) <- sub(paste0("^", colnames(cell.var)[1]), "",
                                   names(scores))
            }
            return(scores)
          })
