###########################    Matrix manipulation    ##########################
#' Sparsen a matrix
#'
#' @description
#' Handy way to convert an object to a sparse dgCMatrix
#'
#' @param mat a dense matrix (or any object that can be converted to a dgCMatrix)
#' @return a sparse dgCMatrix
#' @seealso \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}
#' @keywords internal
#' @noRd

as.dgcmatrix <- function(mat) {
  return(as(mat, "dgCMatrix"))
}

#' Matrix indexing
#'
#' @description
#' Switch between rows & columns indices and vector-like indices of matrices.\cr
#' \cr
#' \code{rowcol2idx} get indices out of rows and columns\cr
#' \code{idx2col} get columns out of indices\cr
#' \code{idx2row} get rows out of indices
#'
#' @param height height of the matrix (number of rows)
#'
#' @return a vector of indices
#' @name matrix-indexing
#' @export

#' @param row.idx indices of rows
#' @param col.idx indices of columns
#' @rdname matrix-indexing
rowcol2idx <- function(rows.idx, cols.idx, height) {
  return((cols.idx - 1) * height + rows.idx)
}

#' @param idx vector-like, matrix-wide indices
#' @rdname matrix-indexing
idx2col <- function(idx, height) {
  return((idx -1) %/% height + 1)
}

#' @rdname matrix-indexing
idx2row <- function(idx, height) {
  return((idx - 1)  %% height + 1)
}

#' Matrix sorting
#'
#' @description
#' Sort a matrix in a row- or column-wise manner, or check that it is\cr
#' \cr
#' \code{rowSort} sort each row independently\cr
#' \code{colSort} sort each column independently\cr\cr
#' \code{rowSorted} check that each row is independently sorted\cr
#' \code{colSorted} check that each column is independently sorted
#'
#' @inheritParams base::sort
#' @param mat matrix to sort or to check
#' @param by matrix to sort by. If \code{NULL}, sort \code{mat} by itself (default)
#' @param ncol number of desired columns for the sorted \code{mat}. If
#' \code{NULL}, keep original dimensions of \code{mat} (default)
#'
#' @return a sorted matrix of \code{ncol} columns
#' @note If \code{mat} and \code{by} have different dimensions, results might be
#' incorrect
#' @name matrix-sorting
#' @export

#' @rdname matrix-sorting
rowSort <- function(mat, by = NULL, ncol = NULL, decreasing = FALSE) {
  by <- by %||% mat
  ncol <- ncol %||% ncol(mat)
  return(.sort.mat(mat, by, byrow = TRUE, ncol = ncol, decreasing = decreasing))
}

#' @rdname matrix-sorting
colSort <- function(mat, by = NULL, ncol = NULL, decreasing = FALSE) {
  by <- by %||% mat
  ncol <- ncol %||% ncol(mat)
  return(.sort.mat(mat, by, byrow = FALSE, ncol = ncol, decreasing = decreasing))
}

#' @keywords internal
.sort.mat <- function(mat, by, byrow, ncol, decreasing) {
  oby <- c(col, row)[[byrow+1]]
  decreasing <- c(FALSE, decreasing)   # c(row or col indices, values)
  return(matrix(mat[order(oby(by), by, decreasing = decreasing, method = "radix")], byrow = byrow, ncol = ncol))
}

#' @rdname matrix-sorting
rowSorted <- function(mat, decreasing = FALSE) {
  o <- 1:ncol(mat)
  if(decreasing) { o <- rev(o) }
  return( .is.sorted.mat(mat = mat[,o,drop=F], MARGIN = 1) )
}
#' @rdname matrix-sorting
colSorted <- function(mat, decreasing = FALSE) {
  o <- 1:nrow(mat)
  if(decreasing) { o <- rev(o) }
  return( .is.sorted.mat(mat = mat[o,,drop=F], MARGIN = 2) )
}

#' @keywords internal
.is.sorted.mat <- function(mat, MARGIN) {
  return( !apply(mat, MARGIN, is.unsorted, strictly = FALSE) )
}


#' Matrix symmetrizing
#'
#' @description
#' Create a sparse matrix using i, j and x (row indices, column indicies and
#' values respectively) then symmetrizes it with low memory footprint\cr
#' \cr
#' \code{symmetrize.pmax.sparse} symmetrize matrix with max(m[i,j], m[j,i])\cr
#' \code{symmetrize.pmin.sparse} symmetrize matrix with min(m[i,j], m[j,i])
#'
#' @inheritParams matrix-indexing
#' @param i row indices (1-based)
#' @param j column indices (1-based)
#' @param x values such that \code{m[i,j] = x}
#'
#' @return a symmetric sparse dgCMatrix of size height x height
#'
#' @seealso \code{\link{SymmetrizeKnn}}
#' @name matrix-symmetrize
#' @export

#' @rdname matrix-symmetrize
symmetrize.pmax.sparse <- function(i, j, x, height) {
  return(.symmetrize.sparse(i, j, x, height, pmax))
}

#' @rdname matrix-symmetrize
symmetrize.pmin.sparse <- function(i, j, x, height) {
  return(.symmetrize.sparse(i, j, x, height, pmin))
}

#' @keywords internal
#' @noRd
.symmetrize.sparse <- function(i, j, x, height, p.min.max) {
  # remove any zero
  i.not.zero <- which(x > 0)
  x <- x[i.not.zero]
  i <- i[i.not.zero]
  j <- j[i.not.zero]

  # keep diagonal values aside
  i.not.diag <- which(i != j)
  i_diag <- i[-i.not.diag]
  j_diag <- j[-i.not.diag]
  x_diag <- x[-i.not.diag]

  i <- i[i.not.diag]
  j <- j[i.not.diag]
  x <- x[i.not.diag]

  # build 2 matrices c(indexes, values) and
  # c(transposed indexes, values)
  skelet.mat <- matrix(c(rowcol2idx(i, j, height), x), ncol = 2)
  skelet.maT <- matrix(c(rowcol2idx(j, i, height), x), ncol = 2)

  # find and add missing indexes
  i.missing <- which(! skelet.maT[,1] %in% skelet.mat[,1])
  skelet.mat <- rbind(skelet.mat, skelet.maT[i.missing,])
  i.missing <- which(! skelet.mat[,1] %in% skelet.maT[,1])
  skelet.maT <- rbind(skelet.maT, skelet.mat[i.missing,])

  # order matrices by index to prepare pmax/pmin
  skelet.mat <- skelet.mat[order(skelet.mat[,1]),]
  skelet.maT <- skelet.maT[order(skelet.maT[,1]),]
  # pmax: keep max values to ensure minimums are dropped
  # not necessary when min values are only zeros
  # pmin: keep min values
  skelet.mat[,2] <- p.min.max(skelet.mat[,2], skelet.maT[,2])

  sym.mat <- sparseMatrix(i=c(idx2row(skelet.mat[,1], height), i_diag),
                          j=c(idx2col(skelet.mat[,1], height), j_diag),
                          x=c(skelet.mat[,2], x_diag),
                          dims = rep(height,2))
  return(sym.mat)

}
################################################################################
#############################    NN symmetrize     #############################
#' Symmetrize a nearest neighbours graph
#'
#' @description
#' The approximate nature of the nearest neighbour search algorithm used to
#' compute the knn graph makes the resulting adjacency matrix asymmetric.
#' Those functions symmetrize a \code{Graph} or a \code{Neighbor} object.
#'
#' @param object a \code{Seurat}, \code{Graph} or \code{Neighbor} object
#' @param graph.name name of a \code{Graph} or \code{Neighbor} instance stored in
#' the \code{Seurat object}.
#' @param use.max by default, use the maximum value in case of discrepancy
#' between m[i,j] and m[j,i]. Set to \code{FALSE} to use the minimum value.
#' @param assay name of the assay to store in the output \code{Graph}
#'
#' @return the Seurat object with a new \code{Graph} instance or a
#' \code{dgCMatrix} representing the \code{Graph} itself
#'
#' @export
#'
#' @seealso The classes \code{\link[SeuratObject:Graph]{Graph}} and
#' \code{\link[SeuratObject:Neighbor]{Neighbor}} and
#' \code{\link{symmetrize.pmax.sparse}}

setGeneric("SymmetrizeKnn",
           function(object, graph.name = "RNA_nn", use.max = TRUE, assay = NULL)
             standardGeneric("SymmetrizeKnn"))

#' @export
#' @rdname SymmetrizeKnn
setMethod("SymmetrizeKnn", "Seurat",
          function(object, graph.name = "RNA_nn", use.max = TRUE, assay = NULL) {
            assay <- assay %||% DefaultAssay(object)
            use.max <- use.max %||% FALSE
            knnmat <- SymmetrizeKnn(object = object[[graph.name]], use.max = use.max)
            if (inherits(object[[graph.name]], "Neighbor")) {
              cells <- slot(object = object[[graph.name]], "cell.names") %||% Cells(object)
            } else {
              cells <- colnames(object[[graph.name]])
            }
            rownames(knnmat) <- colnames(knnmat) <- cells
            new.graph = sprintf("%s_symmetric", graph.name)
            object[[new.graph]] <- as.Graph(knnmat)
            slot(object = object[[new.graph]], name = "assay.used") <- assay
            return(object)
          })

#' @export
#' @rdname SymmetrizeKnn
setMethod("SymmetrizeKnn", "Matrix",
          function(object, use.max = TRUE) {
            i <- slot(object = object, name = "i") + 1
            x <- slot(object = object, name = "x")
            p <- slot(object = object, name = "p")
            j <- findInterval(seq(x)-1,p[-1]) + 1
            symmetrize <- c(symmetrize.pmin.sparse,
                            symmetrize.pmax.sparse)[use.max + 1]
            return(
              symmetrize(i = i, j = j, x = x, height = ncol(object))
            )
          })

#' @export
#' @rdname SymmetrizeKnn
# actually useless since Graph objects get Matrix dispatch
setMethod("SymmetrizeKnn", "Graph", getMethod("SymmetrizeKnn", "Matrix"))

# #' @export
# #' @S3method SymmetrizeKnn Graph
# #' @rdname SymmetrizeKnn
# setMethod("SymmetrizeKnn", "Matrix",
#           getMethod("SymmetrizeKnn", "Graph"))

#' @export
#' @rdname SymmetrizeKnn
setMethod("SymmetrizeKnn", "Neighbor",
          function(object, use.max = TRUE) {
            knn.idx  <- slot(object = object, name = "nn.idx")
            knn.dist <- slot(object = object, name = "nn.dist")
            n <- nrow(knn.idx)
            k <- ncol(knn.idx)
            symmetrize <- c(symmetrize.pmin.sparse,
                            symmetrize.pmax.sparse)[use.max + 1]
            return(
              symmetrize(i = rep(1:n, k), j = as.vector(knn.idx),
                         x = as.vector(knn.dist), height = n)
            )
          })


################################################################################
#############################    NN per batches    #############################
#' Calculate number of nearest neighbours between batches
#'
#' @description
#' Calculate number of nearest neighbours between batches out of a knn graph
#'
#' @param object a \code{Seurat} object
#' @param batch.var name of a column in the \code{Seurat} object's metadata
#' containing batch information
#' @param graph.name name of a \code{Graph} or \code{Neighbor} instance stored in
#' \code{object}. When available, prefer the distance based network to the
#' connectivities graph (especially when computed with the UMAP method).
#' @param count.self whether to include self-to-self vertices in the calculation
#'
#' @return a square count matrix of size number of batches (see \strong{Details}
#' section)
#'
#' @export
#' @details
#' The output matrix will likely not be symmetrical. This is due to the approximate
#' nature of the nearest neighbour search algorithm used to compute the knn
#' graph. It must be read by \strong{row}. For instance, the number of times
#' cells of batch 1 have cells of batch 3 in their nn is matrix[1,3]
#'
#' @seealso The classes \code{\link[SeuratObject:Graph]{Graph}} and
#' \code{\link[SeuratObject:Neighbor]{Neighbor}}

setGeneric("GetNeighborsPerBatch",
           function(object, batch.var, graph.name = "RNA_nn", count.self = TRUE)
             standardGeneric("GetNeighborsPerBatch"),
           signature = c("object", "batch.var"))

#' @keywords internal
setMethod("GetNeighborsPerBatch", c("Seurat", "character"),
  function(object, batch.var, graph.name = "RNA_nn", count.self = TRUE) {
    batch.var <- object[[]][, batch.var, drop = FALSE]
    colnames(batch.var) <- colnames(batch.var) %||% "batch"
    GetNeighborsPerBatch(object = object[[graph.name]], batch.var = batch.var,
                         count.self = count.self)
})

#' @keywords internal
setMethod("GetNeighborsPerBatch", c("Seurat", "data.frame"),
  function(object, batch.var, graph.name = "RNA_nn", count.self = TRUE){
    GetNeighborsPerBatch(object = object[[graph.name]], batch.var = batch.var,
                         count.self = count.self)
})

#' @keywords internal
setMethod("GetNeighborsPerBatch", "Graph",
  function(object, batch.var, count.self) {
    knnmat <- as.dgcmatrix(object > 0)
    GetNeighborsPerBatch(object = knnmat, batch.var = batch.var,
                         count.self = count.self)
})

#' @keywords internal
setMethod("GetNeighborsPerBatch", "matrix", getMethod("GetNeighborsPerBatch", "Graph"))

#' @keywords internal
setMethod("GetNeighborsPerBatch", "Neighbor",
  function(object, batch.var, count.self){
    graph <- as.Graph(object)
    GetNeighborsPerBatch(object = graph, batch.var = batch.var,
                         count.self = count.self)
})

#' @keywords internal
setMethod("GetNeighborsPerBatch", "Matrix",
  function(object, batch.var, count.self) {
    batch.var.nm <- colnames(batch.var)[1]
    formula <- as.formula(sprintf("~ 0 + %s", batch.var.nm))
    binmat <- Matrix::sparse.model.matrix(object = formula, data = batch.var)

    # if binmat has rownames, ensure object matrix dimnames are properly ordered
    # else, assume they are ordered likewise
    order.cells <- rownames(binmat) %||% 1:nrow(binmat)
    object <- object[order.cells, order.cells]
    Matrix::diag(object) <- as.integer(count.self)
    object <- drop0(object)
    binmat <- binmat[order.cells, ]
    nn.per.batch <- t(binmat) %*% (object %*% binmat)

    return(as.matrix(nn.per.batch))
})


#' Calculate proportion of nearest neighbours within batches
#'
#' @description
#' Calculate the proportion of nearest neighbours within batches out of a knn
#' graph
#'
#' @inheritParams GetNeighborsPerBatch
#' @param per.batch whether to keep proportions per batch or to aggregate
#' everything
#'
#' @return a vector of length 1 or number of batches, depending on
#' \code{per.batch} argument value
#'
#' @export
#'
#' @seealso \code{\link{GetPropInterBatch}} and \code{\link{GetNeighborsPerBatch}}

GetPropIntraBatch <- function(object, batch.var, graph.name = "RNA_nn",
                              count.self = TRUE, per.batch = TRUE) {
  nn.per.batch <- GetNeighborsPerBatch(object = object, batch.var = batch.var,
                                       graph.name = graph.name,
                                       count.self = count.self)
  nn.intra.batch <- diag(nn.per.batch)
  nn.total <- rowSums(nn.per.batch)
  if (! per.batch) {
    nn.intra.batch <- sum(nn.intra.batch)
    nn.total <- sum(nn.total)
  }
  return(nn.intra.batch / sum(nn.per.batch))
}


#' Calculate proportion of nearest neighbours between batches
#'
#' @description
#' Calculate the proportion of nearest neighbours between batches out of a knn
#' graph
#'
#' @inheritParams GetPropIntraBatch
#' @inherit GetPropIntraBatch return
#'
#' @export
#'
#' @seealso \code{\link{GetPropIntraBatch}} and \code{\link{GetNeighborsPerBatch}}

GetPropInterBatch <- function(object, batch.var, graph.name = "RNA_nn",
                              count.self = TRUE, per.batch = TRUE) {
  prop.intra.batch <- GetPropIntraBatch(object = object, batch.var = batch.var,
                                        graph.name = graph.name,
                                        count.self = count.self,
                                        per.batch = per.batch)
  return(1 - prop.intra.batch)
}

################################################################################
#################################    Checks    #################################
#' Check whether a knn Graph or Neighbor has a constant k value
#'
#' @keywords internal
#' @noRd
is.kconstant <- function(object) {
  UseMethod(generic = "is.kconstant", object = object)
}
#' @keywords internal
#' @noRd
is.kconstant.Neighbor <- function(object) {
  return(TRUE)
}
#' @keywords internal
#' @noRd
is.kconstant.Matrix <- function(object) {
  return(length(unique(rowSums(object > 0))) == 1)
}

#' Check whether a Graph or Neighbor
#'
#' @keywords internal
#' @noRd
could.be.connectivity <- function(object, check.symmetry = T) {
  UseMethod(generic = "could.be.connectivity", object = object)
}
#' @importFrom SeuratObject as.Graph
#' @importFrom Matrix isSymmetric
#' @keywords internal
#' @noRd
could.be.connectivity.Neighbor <- function(object, check.symmetry = T) {
  cat(check.symmetry, "\n")
  knn.dist <- slot(object = object, name = "nn.dist")
  res <- min(knn.dist) >= 0 & max(knn.dist) <= 1 &
    (! check.symmetry || isSymmetric(as.Graph(object)))
  return(res)
}
#' @keywords internal
#' @noRd
could.be.connectivity.Matrix <- function(object, check.symmetry = T) {
  cat(check.symmetry, "\n")
  res <- min(object) >= 0 & max(object) <= 1 &
    (! check.symmetry || isSymmetric(object))
  return(res)
}
