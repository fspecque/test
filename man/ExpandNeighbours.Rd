% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expand_knn_graph.R
\name{ExpandNeighbours}
\alias{ExpandNeighbours}
\alias{ExpandNeighbours,Seurat-method}
\title{Expand knn graph to increase the number of neighbours}
\usage{
ExpandNeighbours(
  object,
  graph.name = "RNA_nn",
  new.graph.name = NULL,
  graph.type = c("distances", "connectivities"),
  k.target = 90L,
  do.symmetrize = FALSE,
  algo = c("dijkstra", "diffusion"),
  which.dijkstra = c("auto", "igraph", "fast", "slow"),
  dijkstra.ncores = 1L,
  dijkstra.tol = 1L,
  diffusion.iter = 26L,
  assay = NULL,
  verbose = TRUE
)

\S4method{ExpandNeighbours}{Seurat}(
  object,
  graph.name = "RNA_nn",
  new.graph.name = NULL,
  graph.type = c("distances", "connectivities"),
  k.target = 90L,
  do.symmetrize = FALSE,
  algo = c("dijkstra", "diffusion"),
  which.dijkstra = c("auto", "igraph", "fast", "slow"),
  dijkstra.ncores = 1L,
  dijkstra.tol = 1L,
  diffusion.iter = 26L,
  assay = NULL,
  verbose = TRUE
)
}
\arguments{
\item{object}{a \code{Seurat}, \code{Graph} or \code{Neighbor} object}

\item{graph.name}{name of a \code{Graph} or \code{Neighbor} instance stored in
the \code{Seurat object}.}

\item{new.graph.name}{name of the return \code{Graph} to store in the
\code{Seurat object}.}

\item{graph.type}{One of "distances" or "connectivities", indicating the type
of metric stored in the graph object.}

\item{k.target}{number of nearest neighbours to reach}

\item{do.symmetrize}{whether to make the input graph symmetric if necessary.
See \strong{Details} section for further explanations}

\item{algo}{One of "dijkstra" or "diffusion". "diffusion" is suited for
connectivity matrices only}

\item{which.dijkstra}{one of "igraph", "fast" or "slow". "auto" (default)
chooses for you. See \strong{Details} section}

\item{dijkstra.ncores}{number of cores to use for Dijkstra's algorithm.
Ignored when \code{which.dijkstra = "igraph"}}

\item{dijkstra.tol}{number of sequential iterations with identical best
neighbours found to consider that Dijkstra's algorithm should be stopped.
Ignored when \code{which.dijkstra = "igraph"}}

\item{diffusion.iter}{maximum number of iterations to reach \code{k.target}}

\item{assay}{name of the assay to store in the output \code{Graph}}

\item{verbose}{whether to print progress messages}
}
\value{
the Seurat object with a new \code{Graph} instance or a
\code{dgCMatrix} representing the \code{Graph} itself
}
\description{
Expand a knn graph to increase the number of nearest neighbours using
Dijkstra's algorithm, or a diffusion algorithm. Dijkstra's algorithm is used
to prepare for LISI score, while diffusion is suited for preparing to compute
the kBET score. Beware that diffusion is designed to work on connectivity
matrices and is not adequate for distance-based networks.
}
\details{
The approximate nature of the nearest neighbour search algorithm used to
compute the knn graph makes the resulting adjacency matrix asymmetric.
It means that \eqn{cell_i} can be a nearest neighbour of \eqn{cell_j} while
the reverse is not true even though the distance between \eqn{cell_i} and
\eqn{cell_j} is lower than the distance between \eqn{cell_j} and some of its
nearest neighbours.

One can choose to keep the graph as it is and consider it as a directed graph
(\code{do.symmetrize = FALSE}).
The alternative solution is to use all computed distances to extend the knn
graph by making the matrix symmetric. Note that connectivity graphs are
already symmetric, so the argument value should have no effect on the result.
}
