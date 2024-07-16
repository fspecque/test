#' Score a corrected or uncorrected PCA using the Local Inverse Simpson Index
#'
#' @description
#' Compute the Local Inverse Simpson's Index (LISI) to estimate batch mixing or
#' cell type mixing (iLISI and cLISI respectively according to Luecken M.D.
#' \emph{et al.}, 2022).
#'
#' @param object A Seurat object
#' @param batch.var The name of the batch variable (must be in the object metadata)
#' @param reduction The name of the reduction to score
#' @param dims The dimensions to consider. All dimensions are used by default
#' @param use.union Whether to use a union of densities (constructs a fake
#' density with pmax). The alternative is to use all density kernels and sum them
#' @param bw.join A function to pick a joint bandwidth. The default is to use the
#' \code{mean()}. You can also chose \code{median} for instance. If you want to
#' use a distinct bandwidth for each batch, just set \code{bw.join = c} (not recommended)
#' @param bw.kernel Bandwidth selector to use for Gaussian kernel density estimation.
#' One of 'nrd0', 'nrd', 'ucv','bcv' or 'sj' (see \link[stats:bw.nrd]{bandwidths}).
#' nrd0 is the default
#' @param ... Additional parameters ('assay', 'search' and 'scale.layer') to pass
#' to Seurat to automatically construct the \code{batch.var} when not provided.
#'
#' @return A single float corresponding to the score of the given reduction
#'
#' @importFrom stats bw.nrd0 bw.nrd bw.ucv bw.bcv bw.SJ density
#' @importFrom dplyr `%>%` group_by summarise across all_of ungroup mutate pick
#'
#' @export
#' @details
#'  The score is computed as follow :
#' \deqn{\sum_{i=1}^{p} \left ( \cfrac{ \bigcup_{b=1}^n D_{bi} }{ Ref_i } * V_i \right )}
#'
#' For a PCA with p dimensions and n batchs. \eqn{D_{bi}} is the vector of density
#' of batch b on the dimension i. \eqn{V_i} is the proportion of variance
#' explained by the \eqn{PC_i}.
#'
#' with \eqn{Ref = \bigcap_{i=1}^n D_{bi}} or \eqn{Ref = \sum_{i=1}^n D_{bi}}
#' when \code{use.union} is \code{TRUE} and \code{FALSE} respectively, for a
#' given dimension \eqn{i}.
#'
#' \code{use.union = FALSE} tends to result in lower score since the
#' reference is the sum of kernel densities instead of maximum
#'
#' @note This score is an adaptation of the principal component regression (PCR)
#' score from Luecken M.D. \emph{et al.}, 2022.
#'
#' @examples
#' \dontrun{
#' obj <- SeuratData::LoadData("pbmcsca")
#' obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj)
#'
#' score.union <- ScoreDensityPC(obj, "Method", "pca", dim = 1:30)
#' score.sum <- ScoreDensityPC(obj, "Method", "pca", dim = 1:30, use.union = FALSE)
#'
#' score.union   # ~ 0.1511
#' score.sum     # ~ 0.0319
#' }
#'
#' @references Luecken, M. D., Büttner, M., Chaichoompu, K., Danese, A.,
#' Interlandi, M., Mueller, M. F., Strobl, D. C., Zappia, L., Dugas, M.,
#' Colomé-Tatché, M. & Theis, F. J. Benchmarking atlas-level data integration in
#' single-cell genomics. Nat Methods 19, 41–50 (2021).
#' \href{https://doi.org/10.1038/s41592-021-01336-8}{DOI}
#'
#' @seealso \code{\link[Seurat]{Layers}} and \code{\link[Seurat]{DefaultAssay}}
#' for \code{...} arguments and \code{\link{ScoreRegressPC}} to compute the PCR
#' score

AddLISIScore <- function(object, batch.var = NULL, cell.var = NULL, reduction = "pca",
                         dims = NULL, perplexity = 30, tol = 1e-5, do.scale = TRUE,
                         nn.method = c("annoy", "rann", "hnsw"),
                         n.trees = 50,
                         annoy.metric = "euclidean",
                         graph.name = paste0("LISI_perp.", as.integer(perplexity)),
                         nn.eps = 0,
                         l2.norm = FALSE,
                         cache.index = FALSE,
                         verbose = TRUE) {
  i <- c(batch.var %iff% T %||% F, cell.var %iff% T %||% F)
  v <- c(batch.var, cell.var)
  new.names <- paste0(c("LISIbatch_", "LISIcell_"), v, "_", reduction)[i]

  lisi.out <- ScoreLISI(object = object, batch.var = batch.var,
                        cell.var = cell.var, reduction = reduction,
                        dims = dims, perplexity = perplexity, tol = tol,
                        nn.method = nn.method, n.trees = n.trees,
                        annoy.metric = annoy.metric, nn.eps = nn.eps,
                        l2.norm = l2.norm, cache.index = cache.index,
                        verbose = verbose, return.graph = TRUE)

  new.names <- setNames(v, new.names)
  new.names <- new.names[new.names %in% colnames(lisi.out$lisi)]

  object <- AddMetaData(object = object,
                        metadata = lisi.out$lisi %>% rename(all_of(new.names)))
  v <- unname(new.names)
  graph.name %iff% {object[[graph.name]] <- lisi.out$graph}
  lisi.score <- lisi.out$lisi %>% summarise(across({{ v }}, median))

  n <- object[[]] %>% summarize(across({{ v }}, n_distinct)) %>% unlist()
  if (do.scale) {
    lisi.score[batch.var] <- (lisi.score[batch.var] - 1) / (n[batch.var] - 1)
    lisi.score[cell.var] <- (n[cell.var] - lisi.score[cell.var]) / (n[cell.var] - 1)
  }
  slot <- sprintf("LISI_perp.%s_%s_%s", as.integer(perplexity),
                  reduction %||% "pca", c("unscaled", "scaled")[do.scale + 1])
  Misc(object = object, slot = slot) <- lisi.score %>% rename(all_of(new.names))
  object
}


ScoreLISI <- function(object, batch.var = NULL, cell.var = NULL, reduction = "pca",
                      dims = NULL, perplexity = 30, tol = 1e-5,
                      nn.method = c("annoy", "rann", "hnsw"),
                      n.trees = 50,
                      annoy.metric = "euclidean",
                      nn.eps = 0,
                      l2.norm = FALSE,
                      cache.index = FALSE,
                      verbose = TRUE, ...) {
  reduction <- reduction %||% "pca"
  return.graph <- list(...)[["return.graph"]] %||% FALSE
  has.vars <- c("batch" = !is.logical(batch.var %||% FALSE),
                "cell"  = !is.logical(cell.var  %||% FALSE))
  if (!any(has.vars)) {
    abort(message = sprintf("please specify at lease one of %s or %s",
                            sQuote("batch.var"), sQuote("cell.var")))
  }
  vars <- c("batch" = batch.var, "cell" = cell.var)
  nn.method <- tolower(nn.method)
  nn.method <- match.arg(nn.method)
  graph.name <- graph.name[1] %||% paste0(DefaultAssay(object), "_nn_LISI")

  if (! reduction %in% Reductions(object)) {
    abort(message = sprintf("%s reduction not in object", sQuote(reduction)))
  }
  dimred <- Embeddings(object = object, reduction = reduction)[, dims, drop = FALSE]
  if (! all(dims %in% 1:ncol(dimred))) {
    abort(message = "some dims are out of range")
  }

  mtdt <- object[[]]
  found.vars <- intersect(vars[has.vars[names(vars)]], colnames(mtdt))
  if (length(found.vars) == 0) {
    abort(message = sprintf("%s not found in the colnames of the metadata",
                            paste(sQuote(vars[has.vars[names(vars)]]),
                                  collapse = " and ")))
  }
  if (length(found.vars) < sum(has.vars)) {
    warning(sprintf("%s was not found in the colnames of the metadata, skipping",
                    sQuote(setdiff(vars[has.vars[names(vars)]], found.vars))),
            call. = FALSE, immediate. = TRUE)
  }
  mtdt <- mtdt %>% mutate(across({{ found.vars }}, ~ as.integer(as.factor(.x))))

  g <- FindNeighbors(dimred, return.neighbor = TRUE, k.param = as.integer(3*perplexity),
                     nn.method = nn.method, n.trees = n.trees,
                     annoy.metric = annoy.metric, nn.eps = nn.eps,
                     verbose = verbose, #graph.name = graph.name,
                     l2.norm = l2.norm, cache.index = cache.index)

  nn.idx  <- t(g@nn.idx[,-1])
  nn.dist <- t(g@nn.dist[,-1])
  lisi.installed <- requireNamespace("lisi", quietly = TRUE)
  if (! lisi.installed) {
    if (pb$print.lisi.msg) {
      message("To benefit from the original and much faster implementation of ",
              "the Local Inverse Simpson’s Index (LISI), you can install the ",
              "lisi package (by the same team that developed harmony) using either:",
              "\n\nremotes::install_github('immunogenomics/lisi')",
              "\n# or",
              "\ndevtools::install_github('immunogenomics/lisi')",
              "\n\nThis message will be shown once per session")
      pb$print.lisi.msg <- FALSE
    }
    message("Computing LISI scores..."[verbose], appendLF = F)
    simpson.indexes <- lapply(mtdt[, found.vars, drop=FALSE], .compute.lsi,
                              nn.dist = nn.dist, nn.idx = nn.idx,
                              perplexity = perplexity, tol = tol)
    message("done.\n"[verbose], appendLF = F)
  } else {
    get.lsi <- lisi::compute_simpson_index
    mtdt <- mtdt %>% mutate(across({{ found.vars }}, ~ .x - 1))
    nn.idx <- nn.idx - 1
    message("Computing LISI scores with lisi package..."[verbose], appendLF = F)
    simpson.indexes <- lapply(mtdt[, found.vars, drop=FALSE], function(v) {
      lisi::compute_simpson_index(D = nn.dist, knn_idx = nn.idx,
                                  batch_labels = v, n_batches = length(unique(v)),
                                  perplexity = perplexity, tol = tol)
    })
    message("done.\n"[verbose], appendLF = F)
  }
  simpson.indexes <- 1 / as.data.frame(simpson.indexes) # inverse SI
  rownames(simpson.indexes) <- g@cell.names
  if (return.graph) {
    simpson.indexes <- list(graph = g, lisi = simpson.indexes)
  }
  return(simpson.indexes)
}

.compute.lsi <- function(nn.dist, nn.idx, batch.var, perplexity = 30, tol = 1e-5, ...) {
  ncells = ncol(nn.dist)
  nbatch = length(unique(batch.var))
  proba = rep(0, nrow(nn.dist))
  simpson = rep(0, ncells)
  logU = log(perplexity)

  for (i in 1:ncells) {
    beta <- 1
    betamin <- -Inf
    betamax <- Inf
    cell.nn.dist <- nn.dist[, i, drop=TRUE]
    e <- environment()
    list2env(.Hbeta(cell.nn.dist, beta, proba), e)
    Hdiff <- H - logU
    tries <- 0
    # Step 1: get neighbour probabilities
    while(abs(Hdiff) > tol && tries < 50) {
      if (Hdiff > 0){
        betamin <- beta
        if (is.infinite(betamax)) {beta <- beta * 2}
        else {beta <- (beta + betamax) / 2}
      } else{
        betamax <- beta
        if (is.infinite(betamin)) {beta <- beta /2}
        else {beta <- (beta + betamin) / 2}
      }

      list2env(.Hbeta(cell.nn.dist, beta, proba), e)
      Hdiff <- H - logU
      tries <- tries + 1
    }

    if (H == 0) {
      simpson[i] <- -1
    }

    # Step 2: compute Simpson's Index
    for (b in 1:nbatch) {
      q <- which(batch.var[nn.idx[,i]] == b)  # indices of cells belonging to batch (b)
      if (length(q) > 0) {
        sumP <- sum(proba[q])
        simpson[i] <- simpson[i] + sumP^2
      }
    }
  }
  return(simpson)
}

.Hbeta <- function(cell.nn.dist, beta, proba) {
  proba <- exp(- cell.nn.dist * beta)
  sumP <- sum(proba)
  if (sumP == 0) {
    H <- 0
    proba <- rep(0, length(cell.nn.dist))
  } else {
    H <- log(sumP) + beta * sum(cell.nn.dist * proba) / sumP
    proba <- proba / sumP
  }
  return(list("H" = H, "proba" = proba))
}
