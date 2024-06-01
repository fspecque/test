#' Score a corrected or uncorrected PCA to estimate batch mixing
#'
#' @description
#' A batch-wise kernel density estimate is computed for each dimension of the PCA.
#' Then, the global area shared by all kernels is calculated and divided by a
#' reference area (see \strong{Arguments} and \strong{Details}) as a proxy to
#' batch mixing. Finally, the aforementioned proportions are weighted by each
#' dimension's contribution to variance.
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

ScoreDensityPC <- function(object, batch.var=NULL, reduction = "pca", dims=NULL,
                           use.union = TRUE, bw.join = mean,
                           bw.kernel = c('nrd0', 'nrd', 'ucv','bcv', 'sj'), ...) {
  bw.kernel <- tolower(bw.kernel)
  bw.kernel <- match.arg(bw.kernel)
  bw.kernel <- switch (bw.kernel, 'nrd0' = bw.nrd0, 'nrd' = bw.nrd,
                       'ucv' = bw.ucv, 'bcv' = bw.bcv, 'SJ' = bw.SJ)
  bw.join <- bw.join %||% .Primitive("c")
  if (! is.function(bw.join)) {
    abort(message = sprintf("%s is not a function", sQuote("bw.join")))
  }

  prep.list <- .prep_ScorePC(object = object, batch.var = batch.var,
                             reduction = reduction, dims = dims, ...)
  list2env(x = prep.list, envir = environment())

  # compute joint bandwidth
  dimred_cols <- colnames(dimred)
  joint.bws <- df.score %>% group_by(pick({{ batch.var }})) %>%
    summarise(across({{ dimred_cols }}, bw.kernel)) %>%
    ungroup() %>% mutate(across({{ dimred_cols }}, bw.join)) %>%
    as.data.frame()

  # compute densities per batch and collect `density()$y`
  dens.dfs <- sapply(colnames(dimred), function(col.nm) {
    dens.mm <- range(df.score[, col.nm, drop = TRUE], na.rm = TRUE)
    sapply(unique(df.score[,batch.var]), function(batch.nm) {
      density(df.score[df.score[,batch.var] == batch.nm, col.nm],
              bw = joint.bws[batch.nm, col.nm],
              from = dens.mm[1], to = dens.mm[2])$y
    }, simplify = "array") %>% as.data.frame()
  }, simplify = F)

  # compute proportion of common (~intersect) area
  prop.area <- sapply(dens.dfs, function(dens) {
    inter <- do.call(pmin, dens)
    union <- do.call(pmax, dens)
    if (! use.union) {
      # ue sum instead
      union <- rowSums(dens)
    }
    sum(inter) / sum(union)
  }, simplify = "array")

  return(sum(prop.area * proportions(dimvar)))
}


#' Score a corrected or uncorrected PCA to estimate batch mixing
#'
#' @description
#' Linearly regresses principal components to predict batches as a proxy to
#' batch mixing. The resulting R2 are then weighted by each dimension's
#' contribution to variance.
#'
#' @inheritParams ScoreDensityPC
#' @param all.at.once Whether to regress all the dimensions at once instead of
#' regressing them separately. Note that enabling this behaviour doesn't allow
#' to weight R2s by the contribution to variance
#' @param adj.r2 Whether to use the adjusted R2 instead of the raw R2
#'
#' @inherit ScoreDensityPC return
#' @importFrom stats lm summary.lm
#' @importFrom dplyr %>% group_by summarise across all_of ungroup mutate pick
#'
#' @export
#' @details The linear regression is
#' \deqn{Batch = PC_i} or
#' \deqn{Batch = \sum_{i=1}^{p}PC_i} when \code{all.at.once = TRUE}
#'
#'
#' The score is computed as follow :
#' \deqn{\sum_{i=1}^{p} \left ( R^2_i * V_i \right )} or just
#' \deqn{R^2} when \code{all.at.once = TRUE}
#'
#' For a PCA with p dimensions, \eqn{PC_i} is the principal component i,
#' \eqn{R^2_i} is the R squared coefficient of the linear regression for the
#' dimension i. \eqn{V_i} is the proportion of variance explained by the
#' \eqn{PC_i}.
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
#' score.each <- ScoreRegressPC(obj, "Method", "pca", dim = 1:30)
#' score.all <- ScoreRegressPC(obj, "Method", "pca", dim = 1:30, all.at.once = TRUE)
#'
#' score.each    # ~ 0.0257
#' score.all     # ~ 0.7553
#' }
#'
#' @inherit ScoreDensityPC note
#' @inherit ScoreDensityPC references
#'
#' @seealso \code{\link[Seurat]{Layers}} and \code{\link[Seurat]{DefaultAssay}}
#' for \code{...} arguments, \code{\link{ScoreDensityPC}} for an alternative and
#' \code{\link{ScoreRegressPC.CellCycle}} to regresses PCs by cell cycle scores.

# Note: regressing on batch (PC_ ~ batch) or on reduction (batch ~ PC_) results
#       in same (adj.)r2
ScoreRegressPC <- function(object, batch.var=NULL, reduction = "pca", dims=NULL,
                           # regress.what = c('batch.var', 'reduction'),
                           all.at.once = FALSE, adj.r2 = FALSE, ...) {

  all.at.once <- all.at.once %||% FALSE
  adj.r2 <- adj.r2 %||% FALSE
  prep.list <- .prep_ScorePC(object = object, batch.var = batch.var,
                             reduction = reduction, dims = dims, ...)
  # names(prep.list)[length(prep.list)] <- "df.lm"
  list2env(x = prep.list, envir = environment())

  # construct formula(s)
  idx <- 1:length(dims)  # not all.at.once
  if (all.at.once) {
    idx <- rep(x = 1, times = length(dims))
  }
  formulas <- tapply(X = colnames(dimred), INDEX = idx, FUN = function(col.nm) {
    Y_X <- c(batch.var, paste(col.nm, collapse = "+"))
    # if (regress.what == 'batch.var') { Y_X <- rev(Y_X) }
    as.formula(paste(Y_X, collapse = ' ~ '))
  }, simplify = F)

  # compute linear regression per batch and collect (adj.)r.squared
  regs <- lapply(formulas, FUN = lm, data = df.score)
  r2get <- paste0("adj."[adj.r2], "r.squared")
  r2 <- sapply(lapply(regs, summary.lm), getElement, name = r2get, simplify = "numeric")

  return(sum(r2 * proportions(dimvar)))  # if all.at.once, equivalent to return(r2)
}

#' Score a corrected or uncorrected PCA to estimate the contribution of S and G2M
#' scores to variance
#'
#' @description
#' Linearly regresses S and G2M scores to predict principal components. The
#' resulting R2 are then weighted by each dimension's contribution to variance.
#'
#' @inheritParams ScoreRegressPC
#' @param s.var The name of the S phase score variable (must be in the object metadata)
#' @param g2m.var The name of the G2M phase score variable (must be in the object metadata)
#'
#' @inherit ScoreDensityPC return
#' @importFrom stats lm summary.lm
#'
#' @export
#' @details The linear regression is
#' \deqn{PC_i = S_{score} + G2M_{score}}
#'
#' The score is computed as follow :
#' \deqn{\sum_{i=1}^{p} \left ( R^2_i * V_i \right )} or just
#' \deqn{R^2} when \code{all.at.once = TRUE}
#'
#' For a PCA with p dimensions, \eqn{PC_i} is the principal component i,
#' \eqn{R^2_i} is the R squared coefficient of the linear regression for the
#' dimension i. \eqn{V_i} is the proportion of variance explained by the
#' \eqn{PC_i}.
#'
#' @examples
#' \dontrun{
#' obj <- SeuratData::LoadData("pbmcsca")
#' obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj)
#' cc <- CellCycleScoring(JoinLayers(so), s.features = cc.genes.updated.2019$s.genes,
#'                        g2m.features = cc.genes.updated.2019$g2m.genes)[[]]
#' so <- AddMetaData(so, cc[,c("S.Score", "G2M.Score", "Phase")])
#'
#' score.cc <- ScoreRegressPC.CellCycle(obj, "Method", "pca", dim = 1:30)
#'
#' score.cc    # ~ 0.0249
#' }
#'
#' @inherit ScoreDensityPC note
#' @inherit ScoreDensityPC references
#'
#' @seealso \code{\link[Seurat]{Layers}} and \code{\link[Seurat]{DefaultAssay}}
#' for \code{...} arguments and \code{\link{ScoreRegressPC}} to regresses PCs by
#' batch.

ScoreRegressPC.CellCycle <- function(object, batch.var = NULL, reduction = "pca",
                                     dims = NULL, s.var = "S.Score",
                                     g2m.var = "G2M.Score", adj.r2 = FALSE, ...) {
  s.var <- s.var %||% "S.Score"
  g2m.var <- g2m.var %||% "G2M.Score"
  prep.list <- .prep_ScorePC(object = object, batch.var = batch.var,
                             s.var = s.var, g2m.var = g2m.var,
                             reduction = reduction, dims = dims, ...)
  list2env(x = prep.list, envir = environment())

  # construct formula(s)
  idx <- 1:length(dims)
  covar <- c(s.var, g2m.var)
  formulas <- tapply(X = colnames(dimred), INDEX = idx,
                     FUN = paste, "~", paste(covar, collapse = " + "),
                     simplify = F)

  # compute linear regression per batch and collect (adj.)r.squared
  regs <- lapply(lapply(formulas, FUN = as.formula), lm, data = df.score)

  r2get <- paste0("adj."[adj.r2], "r.squared")
  r2 <- sapply(lapply(regs, summary.lm), getElement, name = r2get, simplify = "numeric")

  return(sum(r2 * proportions(dimvar)))
}

#' @importFrom SeuratObject Reductions Embeddings DefaultAssay Layers
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
.prep_ScorePC <- function(object, batch.var = NULL, reduction = "pca",
                          dims = NULL, s.var = NULL, g2m.var = NULL, ...) {
  reduction <- reduction %||% "pca"
  varargs <- list(...)
  if (! reduction %in% Reductions(object)) {
    abort(message = sprintf("%s reduction not in object", sQuote(reduction)))
  }
  dimred <- Embeddings(object = object, reduction = reduction)
  dims <- unique(dims %||% 1:ncol(dimred))
  if (! all(dims %in% 1:ncol(dimred))) {
    abort(message = "some dims are out of range")
  }
  dimred <- dimred[, dims, drop=FALSE]
  dimvar <- Reductions(object = object, slot = reduction)@stdev[dims]^2

  idcol <- "cellid"
  df.score <- object[[]] %>% rownames_to_column(var = idcol)
  if ((s.var %iff% ! all(s.var %in% colnames(df.score))) %||% FALSE) {
    rlang::abort(
      message = sprintf("%s with S phase score not in colnames of metadata",
                        sQuote(s.var)))
  }
  if ((g2m.var %iff% ! all(g2m.var %in% colnames(df.score))) %||% FALSE) {
    abort(
      message = sprintf("%s with G2M phase score not in colnames of metadata",
                        sQuote(g2m.var)))
  }
  batch.var <- batch.var %||% {
    assay <- varargs[["assay"]] %||% DefaultAssay(object) %||% "RNA"
    search <- varargs[["search"]] %||% varargs[["layer"]] %||%
      varargs[["layers"]] %||% "count"
    scale.layer <- varargs[["scale.layer"]] %||% "scale.data"
    df.score <- df.score %>% left_join(
      Seurat:::CreateIntegrationGroups(object = object[[assay]],
                                       layers = Layers(object, search = search),
                                       scale.layer = scale.layer) %>%
        rownames_to_column(idcol), by = idcol
    )
    "group"
  }
  if (! batch.var %in% colnames(df.score)) {
    abort(message = sprintf("%s not in colnames of metadata", sQuote(batch.var)))
  }
  df.score[,batch.var] <- as.numeric(as.factor(df.score[,batch.var]))
  df.score <- df.score %>% left_join(
    as.data.frame(dimred) %>% rownames_to_column(var = idcol),
    by = idcol
  )
  return(list(batch.var = batch.var, dimred = dimred, dimvar = dimvar,
              dims = dims, df.score = df.score))
}
