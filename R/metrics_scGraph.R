#' Score a dimension reduction by correlating normalised centroid distances
#'
#' @description
#' Compute a score based on correlation between batch-wise and dataset-wise
#' normalised centroid distances between cell types. A PCA is computed for each
#' batch and is compared to the reduction output by an integration method.
#'
#' @param object A Seurat object
#' @param cell.var The name(s) of the column(s) with cell type label variable
#' (must be in the object metadata). Multiple column names are accepted
#' @param batch.var The name of the batch variable (must be in the object metadata).
#' Can be omitted if the Seurat object contains multiple layers
#' @param reduction The name of the reduction(s) to score. Multiple names are
#' allowed
#' @param ndims.use Number of dimensions from \code{reduction} to compute
#' centroids with. If the lengths of \code{reduction} and \code{ndims.use} do
#' not match, \code{ndims.use} is recycled or trimmed to have the same length
#' as \code{reduction}.
#' @param ndims.batch Number of principal components to compute for batch-wise
#' PCAs (and to compute centroids with)
#' @param nfeatures Number of features to scale before computing batch-wise PCAs
#' @param approx whether to use truncated singular value decomposition to
#' approximate batch-wise PCAs.
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from each
#' end of PCs or reductions when computing centroids.
#' @param min_cells_type Threshold below which cell types containing too few
#' cells will be excluded
#' @param min_cells_batch Threshold below which batches containing too few
#' cells will be excluded
#' @param cor_method method to compute correlation. One of "weighted_pearson",
#' "pearson" and "spearman".
#' @param assay the name of the assay to use
#' @param verbose whether to print progress messages
#'
#' @return \code{ScoreScGraph}: a named list, with on element per reduction.
#' Each element corresponds to a named vector of raw scores (names are cell type
#' variables)
#'
#' \code{AddScoreScGraph}: the updated Seurat \code{object} with the scGraph
#' raw score(s) set for the integrations.
#'
#' @importFrom SeuratObject Reductions Layers AddMetaData
#' @importFrom stats setNames
#'
#' @export
#' @details
#' For each batch, a PCA is computed using batch-specific variable features.
#' Then, a centroid (trimmed mean) is computed for each cell type and each PC.
#' Cell-type to cell-type centroid distances are then computed for each PC. For
#' each batch, distances are then normalised (divided by the largest distance)
#' in a cell-type-wise manner and averaged across batches to obtain a single
#' mean distance value between cell types.
#'
#' The procedure is then partially repeated once (computation of centroids and
#' distances and normalisation of such distances) on the (un)integrated
#' dimensional reduction of full dataset.
#'
#' Finally, the correlations between the first and the second distance values
#' are computed for each cell-type. The score is obtained by averaging the
#' correlation values. It is bounded by -1 (worse) and 1 (best).
#'
#' @note \code{\link{ScaleScores}()} scales the scores with
#' \deqn{\displaystyle Score_{scaled} = \frac{Score_{raw} + 1}{2}}.
#' @note This score is based on a preprint (see \strong{References})
#'
#' @references Wang H., Leskovec J., Regev A. Metric Mirages in Cell Embeddings.
#' bioRxiv 2024.04.02.587824 [Preprint] (2024).
#' \href{https://doi.org/10.1101/2024.04.02.587824}{DOI}
#' @rdname score-scgraph

ScoreScGraph <- function(object, cell.var, batch.var = NULL,
                         reduction = Reductions(object), ndims.use = NULL,
                         ndims.batch = 10L, nfeatures = 1e3L,
                         approx = FALSE, trim = .2, min_cells_type = 20L,
                         min_cells_batch = 20L,
                         cor_method = c("weighted_pearson", "pearson", "spearman"),
                         assay = NULL, verbose = TRUE) {

  assay <- assay %||% DefaultAssay(object)
  reduction <- reduction %||% Reductions(object)
  cor_method <- tolower(cor_method)
  cor_method <- match.arg(cor_method)

  reduction.in <- intersect(reduction, Reductions(object))
  if (length(reduction.in) == 0) {
    abort('Cannot find any of the provided DimReduc objects')
  }
  if (length(missing_reduc <- setdiff(reduction, reduction.in)) > 0) {
    warning(sprintf('Ignoring unknown DimReduc object(s) %s',
                    paste0(sQuote(missing_reduc), collapse = ", ")),
            call. = F, immediate. = T)
  }
  max_dims <- sapply(reduction.in, function(reduc)
    ncol(Reductions(object, slot = reduc)), simplify = TRUE, USE.NAMES = TRUE)
  ndims.use <- ndims.use %iff%
    setNames(rep(ndims.use, length = length(reduction)), reduction) %||%
    max_dims
  ndims.use <- ndims.use[reduction.in]
  if (any(ndims.use > max_dims)) {
    for (i in which(ndims.use > max_dims)) {
      warning(sprintf("%s has %d dimensions, setting `ndims.use` to %d (instead of %d)",
                      sQuote(reduction.in[i]), max_dims[[i]], max_dims[[i]], ndims.use[[i]]),
              call. = FALSE, immediate. = TRUE)
      ndims.use[[i]] <- max_dims[[i]]
    }
  }

  batch.var <- batch.var[1] %||% {
    group <- CreateIntegrationGroups(
      object[[assay]],
      layers = Layers(object[[assay]], search = 'data'),
      scale.layer = 'scale.data')
    object <- AddMetaData(object, group)
    'group'
  }
  # return `batch.var` and `batches_keep`
  .scgraph_filter_batches(object = object, batch.var = batch.var,
                          min_cells_batch = min_cells_batch, assay = assay)
  # return `cell.var` and `cells_keep`
  .scgraph_filter_cells(object = object, cell.var = cell.var,
                        min_cells_type = min_cells_type)
  # case not covered: less that 2 cell types per batch, causes error when
  #                   computing distances
  # -> covered right before computing distances

  scores <- .scGraph(object = object, cell.var = cell.var, batch.var = batch.var,
                     cells_keep = cells_keep, batches_keep = batches_keep,
                     reduction = reduction.in, ndims.use = ndims.use,
                     ndims.batch = ndims.batch, nfeatures = nfeatures,
                     approx = approx, trim = trim, min_cells_type = min_cells_type,
                     min_cells_batch = min_cells_batch,
                     cor_method = cor_method, assay = assay, verbose = verbose)

}

#' @param integration name of the integration(s) to score. Should match the
#' length of \code{reduction} argument.
#' @export
#' @rdname score-scgraph
AddScoreScGraph <- function(object, integration,
                            cell.var, batch.var = NULL,
                            reduction, ndims.use = NULL,
                            ndims.batch = 10L, nfeatures = 1e3L,
                            approx = FALSE, trim = .2, min_cells_type = 20L,
                            min_cells_batch = 20L,
                            cor_method = c("weighted_pearson", "pearson", "spearman"),
                            assay = NULL, verbose = TRUE) {
  if (length(integration) != length(reduction)) {
    abort("Please provide as many integration names as there are reductions to score.")
  }
  if (length(integration) == 0) {
    abort("Neither `integration` nor `reduction` can be empty")
  }

  integration <- setNames(integration, reduction)

  scores <- ScoreScGraph(object = object, cell.var = cell.var,
                         batch.var = batch.var, reduction = reduction,
                         ndims.use = ndims.use, ndims.batch = ndims.batch,
                         nfeatures = nfeatures, approx = approx, trim = trim,
                         min_cells_type = min_cells_type, assay = assay,
                         min_cells_batch = min_cells_batch,
                         cor_method = cor_method, verbose = verbose)

  object <- check_misc(object)
  for (i in 1:length(scores)) {
    score.names <- paste("scGraph", names(scores[[i]]), sep = '_')
    reduc.name <- names(scores)[i]
    for (j in 1:length(score.names)) {
      object <- SetMiscScore(object, integration = integration[[reduc.name]],
                             score.name = score.names[j],
                             score.value = scores[[i]][[j]],
                             class = "numeric")
    }
  }
  return(object)

}

#' Workhorse of scGraph score computation
#'
#' @importFrom SeuratObject Reductions JoinLayers Embeddings Key
#' @importFrom pbapply pblapply pbsapply pbmapply
#' @importFrom Seurat SplitObject FindVariableFeatures ScaleData RunPCA
#' @importFrom purrr imap map2
#' @importFrom matrixStats rowMaxs
#' @importFrom dplyr %>% bind_rows group_by summarise across everything
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom stats cor
#' @keywords internal
#' @noRd
.scGraph <- function(object, cell.var, cells_keep, batch.var = NULL,
                     batches_keep, reduction = Reductions(object),
                     ndims.use = NULL, ndims.batch = 10L, nfeatures = 1e3L,
                     approx = FALSE, trim = .2, min_cells_type = 20L,
                     min_cells_batch = 20L,
                     cor_method = c("weighted_pearson", "pearson", "spearman"),
                     assay = NULL, verbose = TRUE) {
  lapply_ <- ifelse(verbose, pblapply, lapply)
  sapply_ <- ifelse(verbose, pbsapply, sapply)
  mapply_ <- ifelse(verbose, pbmapply, mapply)

  if (is_installed("distances")) {
    dist_fun <- function(x, y) as.matrix(distances::distances(x, id_variable = y))
  } else {
    dist_fun <- function(x, y) as.matrix(stats::dist(column_to_rownames(x, y)))
  }

  message("Scaling and computing PCA for each batch\n"[verbose], appendLF = F)
  batch_pca <- lapply_(
    SplitObject(JoinLayers(object, assay = assay, layers = "data"), split.by = batch.var)[batches_keep], function(batch_obj) {
      suppressWarnings({
        batch_obj <- FindVariableFeatures(batch_obj, nfeatures = nfeatures, verbose = FALSE)
        batch_obj <- ScaleData(batch_obj, verbose = FALSE)
        batch_obj <- RunPCA(batch_obj, npcs = ndims.batch,
                            reduction.name = "pcascgraph",
                            reduction.key = "pcascgraph_", approx = approx,
                            verbose = FALSE)
      })
      return(Embeddings(batch_obj, reduction = "pcascgraph"))
    })

  message("Computing centroids for each batch and cell type\n"[verbose], appendLF = F)
  centroids_batch <- sapply_(cell.var, function(cell.var_) {
    cell.var.col <- object[[]][, cell.var_, drop = FALSE]
    lapply(batch_pca, .trimmed_mean_centroids, reduction.key = "pcascgraph_",
           cell.var.col = cell.var.col, cell_type_keep = cells_keep[[cell.var_]],
           npcs = ndims.batch, trim = trim)
  }, simplify = FALSE)

  for (cell.var_ in cell.var) {
    l <- sapply(centroids_batch[[cell.var_]], nrow)
    too_small <- l %in% c(0, 1)
    if (all(too_small)) {
      warning(sprintf("`cell.var` %s skipped because it contain less than 2 cell types per batch",
                      sQuote(cell.var_)),
              call. = FALSE, immediate. = TRUE)
      cell.var <- setdiff(cell.var, cell.var_)
      centroids_batch <- centroids_batch[cell.var]
      cells_keep <- cells_keep[cell.var]
    } else {
      if (any(too_small)) {
        warning(sprintf("batches %s ignored because they contain less than 2 cell types per batch for `cell.var` %s",
                        paste(sQuote(names(l)[too_small]), collapse = ", "), cell.var_),
                call. = FALSE, immediate. = TRUE)
        centroids_batch[[cell.var_]] <- centroids_batch[[cell.var_]][! too_small]
      }
    }
  }
  if (length(centroids_batch) == 0) {
    abort(paste("Nothing left to score, try to increase `min_cells_type` and/or",
                "`min_cells_batch`. Make sure that each batch contains at",
                "least 2 cell types"))
  }

  dist_centroids_batch <- imap(centroids_batch, map2, dist_fun)
  dist_centroids_batch <- lapply(dist_centroids_batch, lapply, function(x) {
    as.data.frame(sweep(x, MARGIN = 2, rowMaxs(x), FUN = `/`)) %>%
      rownames_to_column("celltype")
  })
  dist_centroids_batch <- lapply(dist_centroids_batch, function(x)
    x %>% bind_rows() %>%
    group_by(celltype) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    column_to_rownames("celltype") %>% as.data.frame()
  )

  dist_centroids_batch <- mapply(function(mat, cells) mat[cells, cells],
                                dist_centroids_batch, cells_keep,
                                SIMPLIFY = FALSE)

  message("Computing centroids for each reduction and cell type\n"[verbose], appendLF = F)
  centroids_integ <- sapply_(reduction, function(reduc) {
    reduc_obj <- Reductions(object, slot = reduc)
    sapply(cell.var, function(cell.var_) {
      cell.var.col <- object[[]][, cell.var_, drop = FALSE]
      .trimmed_mean_centroids(reduction = Embeddings(reduc_obj),
                              reduction.key = Key(reduc_obj),
                              cell.var.col = cell.var.col,
                              cell_type_keep = cells_keep[[cell.var_]],
                              npcs = ndims.use[[reduc]], trim = trim)
    }, simplify = FALSE)
  }, simplify = FALSE)

  dist_centroids_integ <- lapply(centroids_integ, imap, dist_fun)
  dist_centroids_integ <- lapply(dist_centroids_integ, lapply, function(x) {
    as.data.frame(sweep(x, MARGIN = 2, rowMaxs(x), FUN = `/`))
  })

  message(paste0("Computing correlations using ",
                 sub("_", " ", cor_method), "\n")[verbose], appendLF = F)
  if (cor_method == "weighted_pearson") {
    score <- function(dist_batch, dist_integ, weights) {
      s <- mapply(cor_wt_pearson, x = as.list(dist_batch),
                  y = as.list(dist_integ), W = asplit(weights, 2))
      mean(s, na.rm = TRUE)
    }
    scores <- lapply(dist_centroids_integ, function(dist_integ_reduc) {
      raw_weights <- lapply(dist_integ_reduc, function(dist_mat) {
        W <- 1/as.matrix(dist_mat)
        W[is.infinite(W)] <- 0L
        W
      })
      mapply(score, dist_batch = dist_centroids_batch,
             dist_integ = dist_integ_reduc, weights  = raw_weights)
    })
  } else {
    score <- function(dist_batch, dist_integ, MoreArgs) {
      s <- mapply(cor, x = as.list(dist_batch), y = as.list(dist_integ),
                  MoreArgs = MoreArgs)
      mean(s, na.rm = TRUE)
    }
    scores <- lapply(dist_centroids_integ, function(dist_integ_reduc) {
      MoreArgs <- list(method = cor_method, use = "pairwise.complete.obs")
      mapply(score, dist_batch = dist_centroids_batch,
             dist_integ = dist_integ_reduc,
             MoreArgs = list(MoreArgs = MoreArgs))
    })
  }
  return(scores)
}
#' Compute centroids (trimmed means) per dimension and cell type
#'
#' @importFrom dplyr %>% left_join filter group_by summarize across starts_with
#' @importFrom tibble rownames_to_column
#' @keywords internal
#' @noRd
.trimmed_mean_centroids <- function(reduction, reduction.key, cell.var.col,
                                    cell_type_keep, npcs, trim) {
  cell.var <- colnames(cell.var.col)[1]
  centroids <- as.data.frame(reduction[, 1:npcs, drop = FALSE]) %>%
    rownames_to_column("cellid") %>%
    left_join(cell.var.col %>% rownames_to_column("cellid"), by = "cellid") %>%
    filter(!! data_sym(cell.var) %in% !! cell_type_keep) %>%
    group_by(!! sym(cell.var)) %>%
    summarize(across(starts_with(!! reduction.key),
                     ~ mean(.x, trim = !! trim, na.rm = TRUE))) %>%
    as.data.frame()
  return(centroids)
}

#' Calculate weighted pearson correlation value
#' @keywords internal
#' @noRd
cor_wt_pearson <- function(x, y, W) {
  idx <- which(!is.na(x) & !is.na(y))
  if(length(idx) > 0) {
    x <- x[idx]
    y <- y[idx]
    W <- W[idx]
    W <- W/sum(W)
    W_x <- sum(x * W) / sum(W)
    W_y <- sum(y * W) / sum(W)
    covar <- sum(W * (x - W_x) * (y - W_y))
    var_x <- sum(W * (x - W_x) ** 2)
    var_y <- sum(W * (y - W_y) ** 2)
    weighted_pearson_corr <- covar / sqrt(var_x * var_y)
  } else {
    weighted_pearson_corr <- NA
  }
  return(weighted_pearson_corr)
}

#' @importFrom dplyr %>% count filter pull
#' @importFrom tibble rownames_to_column
#' @keywords internal
#' @noRd
.scgraph_filter_batches <- function(object, batch.var, min_cells_batch, assay) {
  if (! batch.var %in% colnames(object[[]])) {
    abort(sprintf("`batch.var` %s not found in the colnames of the metadata", batch.var))
  }
  batches <- unique(object[[]][, batch.var, drop = TRUE])
  batches_keep <- object[[]] %>% count(!! sym(batch.var)) %>%
    filter(n >= min_cells_batch) %>%
    pull(!! sym(batch.var))

  if (length(batches_keep) == 0) {
    abort("None of the batches contained enough cells. Try to reduce `min_cells_batch`")
  }
  if (length(missing_batches <- setdiff(batches, batches_keep)) > 0) {
    warning(sprintf("%s batches skipped because they contain less than %d cells",
                    paste0(sQuote(missing_batches), collapse = ", "),
                    min_cells_batch))
  }
  list2env(list(batch.var = batch.var, batches_keep = batches_keep),
           envir = parent.frame())
}

#' @importFrom dplyr %>% count filter pull
#' @importFrom tibble rownames_to_column
#' @keywords internal
#' @noRd
.scgraph_filter_cells <- function(object, cell.var, min_cells_type) {
  cell.var.in <- cell.var %in% colnames(object[[]])
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
  cells <- lapply(as.list(object[[]][, cell.var, drop = FALSE]), unique)
  cells_keep <- sapply(cell.var, function(col.name) {
    object[[]] %>% count(!! sym(col.name)) %>%
      filter(n >= min_cells_type) %>%
      pull(!! sym(col.name))
  }, USE.NAMES = TRUE, simplify = FALSE)
  cell_types_rm <- lengths(cells_keep) %in% c(0, 1)
  if (all(cell_types_rm)) {
    abort("scGraph requires that at least 2 cell types contain enough cells per `cell.var`. Try to reduce `min_cells_type`")
    # abort("Less than 2 cell types contained enough cells for all `cell.var`. Try to reduce `min_cells_type`")
  }
  if (any(cell_types_rm)) {
    warning(sprintf("`cell.var` %s skipped because they contain less than 2 cell types with more cells than `min_cells_type` = %d",
                    paste0(names(cells_keep)[cell_types_rm], collapse = ", "), min_cells_type),
            call. = FALSE, immediate. = TRUE)
    cell.var <- cell.var[!cell_types_rm]
    cells <- cells[!cell_types_rm]
    cells_keep <- cells_keep[!cell_types_rm]
  }

  if ((n_missing_types <- sum(lengths(cells)) - sum(lengths(cells_keep))) > 0) {
    warning(sprintf("%d cell types skipped because they contain less than %d cells",
                    n_missing_types, min_cells_type),
            call. = FALSE, immediate. = TRUE)
  }
  list2env(list(cell.var = cell.var, cells_keep = cells_keep),
           envir = parent.frame())
}
