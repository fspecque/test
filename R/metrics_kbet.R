# metrics_pca.R -> .prep_MetaDataBatch
#' @include utils.R
#' @include convert_distances.R
#' @include expand_knn_graph.R
#' @include metrics_pca.R
#' @include scores.R
NULL

#' Score an embedding or a knn graph with the kBET test
#'
#' @description
#' Compute scores based on the k-nearest neighbour batch effect test (kBET).
#' Accepts a \code{\link[SeuratObject:DimReduc-class]{DimReduc}} instance or a
#' graph object (\code{\link[SeuratObject:Graph-class]{Graph}} or
#' \code{\link[SeuratObject:Neighbor-class]{Neighbor}}) provided that it
#' contains distances or connectivities.
#'
#' If '\code{what}' is a \code{\link[SeuratObject:DimReduc-class]{DimReduc}}
#' object, a knn graph is computed for \code{k = 50}. Then, (or if '\code{what}'
#' is a \code{\link[SeuratObject:Graph-class]{Graph}} or a
#' \code{\link[SeuratObject:Neighbor-class]{Neighbor}} object), connectivities
#' are derives from 50 nearest neighbours distances. Finally, the kBET test is
#' computed for each cell-type label of all provided \code{cell.var}. kBET
#' measures if local batch label distribution is similar to the global one. The
#' result of kBET is the average test rejection rate (between 0 and 1). The
#' closest to zero, the less bias is attributable to the batch effect.
#'
#' @param object A Seurat object
#' @param batch.var The name of the batch variable (must be in the object
#' metadata)
#' @param cell.var The name of the column with cell type label variable
#' (must be in the object metadata). If \code{NULL} is passed, the kBET score is
#' computed for all cells at once.
#' @param what the name of the dimension reduction or Graph/Neighbor object to
#' score. Must be in the Seurat object.
#' @param graph.type one of 'distances' (default) or 'connectivities', to
#' indicate the type of graph when \code{what} is not a dimension reduction.
#' @param seed.use the value of the seed to obtain reproducible results.
#' \code{NULL} disables to use of a seed
#' @param verbose whether to print progress messages
#' @param assay assay to use. Passed to Seurat to automatically construct the
#' \code{batch.var} when not provided
#' (\code{\link[SeuratObject:DefaultAssay]{DefaultAssay()}} is used by default).
#' Useless otherwise
#' @param layer layer to use. Passed to Seurat to automatically construct the
#' \code{batch.var} when not provided. Useless otherwise
#'
#' @return \code{ScoreKBET}: a list with one element per \code{cell.var}. Each
#' is a named numeric vector of kBET score (floats between 0 and 1), one per
#' cell-label.
#'
#' \code{AddScoreKBET}:  the updated Seurat \code{object} with the mean kBET
#' score(s) set for the integration.
#'
#' @importFrom SeuratObject DefaultAssay as.Neighbor
#' @importFrom Seurat FindNeighbors
#' @importFrom rlang abort !! data_sym
#' @importFrom dplyr %>% group_by pick summarise n n_distinct filter pull
#' @importFrom igraph components graph_from_adjacency_matrix
#' @importFrom Matrix spMatrix
#' @importFrom kBET kBET
#'
#' @export
#' @references Büttner, M., Miao, Z., Wolf, F. A., Teichmann, S. A. & Theis, F.
#' J. A test metric for assessing single-cell RNA-seq batch correction.
#' Nat Methods 16, 43–49 (2018).
#' \href{https://doi.org/10.1038/s41592-018-0254-1}{DOI}
#' @references Luecken, M. D., Büttner, M., Chaichoompu, K., Danese, A.,
#' Interlandi, M., Mueller, M. F., Strobl, D. C., Zappia, L., Dugas, M.,
#' Colomé-Tatché, M. & Theis, F. J. Benchmarking atlas-level data integration in
#' single-cell genomics. Nat Methods 19, 41–50 (2021).
#' \href{https://doi.org/10.1038/s41592-021-01336-8}{DOI}
#' @rdname score-kbet
ScoreKBET <- function(object, batch.var, cell.var, what,
                      graph.type = c("distances", "connectivities"),
                      seed.use = 42L, verbose = TRUE,
                      assay = NULL, layer = NULL) {
  graph.type <- tolower(graph.type)
  graph.type <- match.arg(graph.type)
  cell.var <- if (isFALSE(cell.var)) NULL else cell.var
  assay <- assay %||% DefaultAssay(object)
  idcol <- "cellbarcodeid"
  .prep_MetaDataBatch(object = object, batch.var = batch.var,
                      assay = assay, layer = layer, idcol = idcol)

  if (cell.var %iff% FALSE %||% TRUE) {
    cell.var <- 'identity'
    df.mtdt$identity <- 1
  }
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

  if (inherits(object[[what]], "DimReduc")) {
    object <- FindNeighbors(object, reduction = what, k.param = 50L,
                            return.neighbor = TRUE, compute.SNN = FALSE,
                            graph.name = "knn.kbet", verbose = verbose)
    what <- "knn.kbet"
  }
  if (inherits(object[[what]], "Graph")) {

    object[["knn.kbet"]] <- as.Neighbor(.CutKnn(object[[what]], k.max = 50L,
                                                assay = assay,
                                                verbose = verbose))
    what <- "knn.kbet"
  }
  if (inherits(object[[what]], "Neighbor")) {
    if (graph.type == "connectivities") {
      if (! could.be.connectivity(object[[what]])) {
        abort(sprintf("%s don't seem to be a connectivity matrix", object[[what]]))
      }
    } else {
      object <- GetConnectivities(object, neighbors.use = what,
                                  graph.name = "conn.kbet", assay = assay,
                                  verbose = verbose)
      what <- 'conn.kbet'
    }

  } else {
    abort(sprintf("Need a DimReduc, Graph or Neighbor object, not %s",
                  sQuote(class(object[[what]]))))
  }

  score.list <- sapply(cell.var, function(vec) {
    k0s <- .estimate.k0.kBET(df.mtdt, celltype.var = vec, batch.var = batch.var)
    crosstab <- df.mtdt %>% group_by(pick({{ vec }})) %>%
      summarise(n_cells = n(), n_batches = n_distinct(pick({{ batch.var }})))

    celltypes.keep <- crosstab %>%
      filter(n_cells >= 10, n_batches > 1) %>% pull({{ vec }})

    n.out <- nrow(crosstab) - length(celltypes.keep)
    if(n.out > 0) {
      warning("skipping ", n.out, " cell types (",
              paste(setdiff(crosstab[[vec]], celltypes.keep), collapse = ', '),
              ") because they contain less than 10 cells or a single batch.")
    }

    scores <- sapply(crosstab[[vec]], function(celltype) {
      if (celltype %in% celltypes.keep) {
        cellids <- df.mtdt %>% filter(!!data_sym(vec) == !!celltype) %>%
          pull({{ idcol }})
        cellmask <- rep_len(TRUE, length(cellids))

        graph.obj <- object[[what]][cellids, cellids]
        graph.comp <- components(
          graph_from_adjacency_matrix(graph.obj > 0,
                                      weighted = NULL),
          mode = 'weak')

        k0 <- k0s[[celltype]]
        knn.idx <- matrix(NA_integer_, nrow = length(cellids), ncol = k0)
        if (graph.comp$no > 1) {
          min.comp.size <- 3 * k0
          prop.cells.large.comp <- with(graph.comp, sum(proportions(csize)[csize >= min.comp.size]))
          if (prop.cells.large.comp >= .75) {
            cellmask <- with(graph.comp, membership %in% seq_along(csize)[csize >= min.comp.size])
            graph.obj <- graph.obj[cellmask, cellmask]
          } else {
            return(1)
          }
        }
        graph.obj <- expand.neighbours.diffusion4kBET(conmat = graph.obj, k.min = k0)
        if (get.k(object = graph.obj, FUN = "min") >= k0) {
          knn.idx[cellmask, ] <- get.k.best.neighbours(
            object = graph.obj, k.max = k0, graph.type = "connectivities"
          )
          seed.use %iff% set.seed(seed.use)
          batch.vec <- df.mtdt[match(cellids, df.mtdt[[idcol]]), batch.var]
          score <- kBET(df = spMatrix(nrow = nrow(knn.idx), ncol = k0),
                        batch = batch.vec, k0 = k0, knn = knn.idx,
                        do.pca = FALSE, heuristic = FALSE, plot = FALSE,
                        adapt = FALSE, verbose = verbose)$summary$kBET.observed[1]
          return(score)

        }
      } else {
        return(1)
      }
    }, simplify = "numeric", USE.NAMES = TRUE)
    return(setNames(scores, crosstab[[vec]]))
  }, simplify = FALSE, USE.NAMES = TRUE)

  return(score.list)
}

#' @param integration name of the integration to score
#' @importFrom dplyr %>%
#' @export
#' @rdname score-kbet
AddScoreKBET <- function(object, integration,
                         batch.var, cell.var, what,
                         graph.type = c("distances", "connectivities"),
                         seed.use = 42L, verbose = TRUE,
                         assay = NULL, layer = NULL) {
  scores <- ScoreKBET(object, batch.var = batch.var, cell.var = cell.var,
                      what = what, graph.type = graph.type, seed.use = seed.use,
                      verbose = verbose, assay = assay, layer = layer)


  score.names <- paste("kBET", names(scores), sep = '_')

  scores <- unlist(lapply(scores, mean, na.rm = TRUE), use.names = TRUE)
  object <- check_misc(object)
  for (i in 1:length(scores)) {
    object <- SetMiscScore(object, integration = integration,
                           score.name = score.names[i],
                           score.value = scores[[i]],
                           class = "numeric")
  }
  return(object)
}

#' @importFrom dplyr %>% count pick group_by mutate case_when pull
#' @keywords internal
#' @noRd
.estimate.k0.kBET <- function(df.mtdf, celltype.var, batch.var,
                              max.int = .Machine$integer.max) {
  k0 <- df.mtdf %>%
    count(pick({{ celltype.var }}, {{ batch.var }}) , name = "n_cells") %>%
    group_by(pick({{ celltype.var }})) %>%
    mutate(k0_bound = floor(mean(n_cells) / 4)) %>%
    mutate(k0 = min(max(10, k0_bound), 70)) %>%
    mutate(k0 = case_when(k0 * sum(n_cells) >= max.int ~ floor(max.int / sum(n_cells)),
                          T ~ k0)) %>%
    pull(k0, name = {{ celltype.var }})
}
