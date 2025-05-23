% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics_lisi.R
\name{AddScoreLISI}
\alias{AddScoreLISI}
\alias{ScoreLISI}
\title{Score a dimensionality reduction embedding or knn graph using the Local
Inverse Simpson Index}
\usage{
AddScoreLISI(
  object,
  integration,
  batch.var = NULL,
  cell.var = NULL,
  reduction,
  dims = NULL,
  graph.name,
  graph.type = c("distances", "connectivities"),
  do.symmetrize = TRUE,
  save.graph = TRUE,
  new.graph = NULL,
  perplexity = 30,
  tol = 1e-05,
  do.scale = TRUE,
  largest.component.only = FALSE,
  assay = NULL,
  verbose = TRUE,
  ...
)

ScoreLISI(
  object,
  batch.var = NULL,
  cell.var = NULL,
  reduction,
  dims = NULL,
  graph.name,
  graph.type = c("distances", "connectivities"),
  do.symmetrize = TRUE,
  return.graph = FALSE,
  perplexity = 30,
  tol = 1e-05,
  largest.component.only = FALSE,
  assay = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{A Seurat object}

\item{integration}{name of the integration to score}

\item{batch.var}{The name of the batch variable (must be in the object metadata).
Can be omitted if \code{cell.var} is not \code{NULL}}

\item{cell.var}{The name of the  cell variable (must be in the object metadata).
Can be omitted if \code{batch.var} is not \code{NULL}}

\item{reduction}{The name of the reduction to score. Arguments
\code{reduction} and \code{graph.name} are mutually exclusive}

\item{dims}{The dimensions of \code{reduction} to consider. All dimensions are
used by default. Ignored when scoring a graph}

\item{graph.name}{The name of the knn graph to score. Arguments
\code{reduction} and \code{graph.name} are mutually exclusive}

\item{graph.type}{one of 'distances' or 'connectivities' (not supported yet).
Ignored when scoring a cell embedding}

\item{do.symmetrize}{whether to symmetrize the knn graphs. Set to\code{FALSE}
to disable (not recommended, especially when scoring a knn graph directly.
See \strong{Details})}

\item{save.graph}{whether to save the graph used to compute the LISI score(s)
in the Seurat object}

\item{new.graph}{name of the graph used to compute the LISI score(s).
When \code{new.graph = NULL} (default), a name is constructed depending on
input and arguments. Ignored when \code{save.graph = FALSE}}

\item{perplexity}{third of number of each cell's neighbours. When the value of
perplexity and the number of neighbours (*3) are discrepant, the graph is
adapted. Multiple scores with varying values of perplexity are not comparable,
hence it is recommended to use the same value for each integration to score.}

\item{tol}{Stop the computation of the local Simpson's index when it converges
to this tolerance.}

\item{do.scale}{whether to scale the output LISI values between 0 and 1.}

\item{largest.component.only}{whether to compute the scores on the largest
component or all sufficiently large components (default)}

\item{assay}{the name of the assay to reference in the output Graph object (
when \code{save.graph = TRUE})}

\item{verbose}{whether to print progress messages}

\item{...}{Additional parameters to pass to other methods
(see \strong{Details}).}

\item{return.graph}{whether to return the graph used to compute the LISI score(s)}
}
\value{
\code{ScoreLISI}: a data frame with each cell's raw LISI score, or
a list containing the aforementioned data frame and the graph used to compute
it (\code{return.graph = TRUE}).

\code{AddScoreLISI}: the updated Seurat \code{object}, with cell-wise LISI
scores in the meta data (identical to \code{ScoreLISI}'s output), global
scores in misc and a new Graph object when \code{save.graph = TRUE}.
}
\description{
Compute the Local Inverse Simpson's Index (LISI) to estimate batch mixing or
cell type mixing (iLISI and cLISI respectively according to Luecken M.D.
\emph{et al.}, 2022).

\code{AddScoreLISI} returns an updated Seurat object, while \code{ScoreLISI}
outputs the raw LISI scores for each cell
}
\details{
When scoring a reduction, a knn graph with enough neighbours per cell is
computed. If \code{do.symmetrize = TRUE}, the graph is symmetrized and the k
best neighbours are kept.

When scoring a knn graph, the graph is expanded with Dijkstra's algorithm to
reach enough neighbours per cell. If \code{do.symmetrize = TRUE}, the graph
is symmetrized beforehand. Note that when \code{do.symmetrize = FALSE},
Dijkstra's algorithm is applied on the asymmetric distance matrix, and
components are computed. But each component is then symmetrized and Dijkstra's
algorithm is computed on each of them. Indeed, there is no guarantee that the
cells' number of neighbours is preserved after decomposing the directed graph.
Hence, when \code{do.symmetrize = FALSE} and a graph is input, the graph is
considered as directed only to find components.

In either case, it is recommended to keep \code{do.symmetrize = TRUE}.

For possible additional parameters, see \code{\link[Seurat]{FindNeighbors}}
(when inputting a reduction) or \code{\link{ExpandNeighbours}} (when
inputting a knn graph)
}
\note{
This score is an adaptation of the LISI score as described in Korsunsky
I. \emph{et al.}, 2019 and also used in Luecken M.D. \emph{et al.}, 2022.
}
\references{
Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F.,
Wei, K., Baglaenko, Y., Brenner, M., Loh, P. & Raychaudhuri, S. Fast,
sensitive and accurate integration of single-cell data with Harmony. Nat
Methods 16, 1289–1296 (2019).
\href{https://doi.org/10.1038/s41592-019-0619-0}{DOI}

Luecken, M. D., Büttner, M., Chaichoompu, K., Danese, A.,
Interlandi, M., Mueller, M. F., Strobl, D. C., Zappia, L., Dugas, M.,
Colomé-Tatché, M. & Theis, F. J. Benchmarking atlas-level data integration in
single-cell genomics. Nat Methods 19, 41–50 (2021).
\href{https://doi.org/10.1038/s41592-021-01336-8}{DOI}
}
\seealso{
\code{\link[Seurat]{FindNeighbors}}, \code{\link{ExpandNeighbours}}
and \code{\link{CutKnn}}
}
