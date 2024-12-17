#' Find a clustering that maximises NMI or ARI
#' @description
#' Compute clusters for multiple resolutions and saves in the metadata the
#' clustering result that reaches the maximum NMI and/or ARI value for a given
#' cell-type label variable.
#'
#' @inheritParams Seurat::FindClusters
#' @param object a Seurat object
#' @param graph.name the name of the knn graph to score.
#' @param cell.var The name(s) of the column(s) with cell type label variable
#' (must be in the object metadata). Multiple column names are accepted
#' @param cluster.name a (optionally 'glue') string used as the new metadata
#' column name (see \strong{Details} section)
#' @param resolutions the resolutions to compute clusters for
#' @param optimisation.metric one of "nmi" or "ari" or both (default).
#' The metric(s) to use to check clustering results against the \code{cell.var}.
#'
#' @return the updated seurat object with the new metadata column(s)
#'
#' @details
#' \code{cluster.name} can use the 'glue' syntax to avoid overwriting metadata
#' columns when multiple metrics and/or cell-type label variables are provided.
#' But it can be ignored otherwise. Injectable variables are "graph.name",
#' "cell.var" and "metric" (\strong{not} "\strong{optimisation.metric}"). They
#' must be \strong{flanked by single curly brackets} ("\{" and "\}").
#'
#' For instance, if you prefer to name the clusters with an integration instead
#' of the \code{graph.name}, don't use the glue syntax for it (e.g.
#' "combat_\{cell.var\}_\{metric\}" would work, but any of
#' "\{combat\}_\{cell.var\}_\{metric\}" or
#' "\{integration\}_\{cell.var\}_\{metric\}" would throw an error)
#'
#' @seealso \code{\link[Seurat]{FindClusters}}, \code{\link{ScoreARI}} and
#' \code{\link{ScoreNMI}}
#'
#' @importFrom SeuratObject CheckDots DefaultAssay AddMetaData LogSeuratCommand
#' @importFrom Seurat FindClusters
#' @importFrom purrr imap
#' @importFrom dplyr %>% select rename bind_cols
#' @export
#'
#' @references Luecken, M. D., Büttner, M., Chaichoompu, K., Danese, A.,
#' Interlandi, M., Mueller, M. F., Strobl, D. C., Zappia, L., Dugas, M.,
#' Colomé-Tatché, M. & Theis, F. J. Benchmarking atlas-level data integration in
#' single-cell genomics. Nat Methods 19, 41–50 (2021).
#' \href{https://doi.org/10.1038/s41592-021-01336-8}{DOI}

FindOptimalClusters <- function(object, graph.name = NULL,
                                cell.var = NULL,
                                cluster.name = "{graph.name}_{cell.var}_{metric}",
                                modularity.fxn = 1,
                                initial.membership = NULL,
                                node.sizes = NULL,
                                resolutions = seq(from = .1, to = 2, by = .1),
                                optimisation.metric = c("nmi", "ari"),
                                method = "matrix",
                                algorithm = 1,
                                n.start = 10,
                                n.iter = 10,
                                random.seed = 0,
                                group.singletons = TRUE,
                                temp.file.location = NULL,
                                edge.file.name = NULL,
                                verbose = TRUE, ...) {
  CheckDots(...)
  optimisation.metric <- tolower(optimisation.metric)
  optimisation.metric <- match.arg(optimisation.metric, several.ok = TRUE)
  cell.var %||% abort('Cluster optimisation require at least a cell type label variable')

  graph.name <- graph.name %||% paste0(DefaultAssay(object = object),
                                       "_snn")
  if (!graph.name %in% names(x = object)) {
    stop("Provided graph.name not present in Seurat object")
  }
  if (!is(object = object[[graph.name]], class2 = "Graph")) {
    stop("Provided graph.name does not correspond to a graph object.")
  }

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

  clustering.results <- FindClusters(object = object[[graph.name]],
                                     modularity.fxn = modularity.fxn,
                                     initial.membership = initial.membership,
                                     node.sizes = node.sizes,
                                     resolution = resolutions, method = method,
                                     algorithm = algorithm, n.start = n.start,
                                     n.iter = n.iter, random.seed = random.seed,
                                     group.singletons = group.singletons,
                                     temp.file.location = temp.file.location,
                                     edge.file.name = edge.file.name,
                                     verbose = FALSE, ...)

  ### maximize metric
  cells <- intersect(rownames(object[[]]), rownames(clustering.results))
  clustering.results <- clustering.results[cells,]
  cell.var.vecs <- asplit(object[[]][cells, cell.var, drop = FALSE], MARGIN = 2)

  best.clusters <- imap(cell.var.vecs, function(cell.var.vec, cell.var) {
    message(sprintf("%s\n", cell.var)[verbose], appendLF = FALSE)
    sapply(optimisation.metric, function(metric) {
      FUN <- switch (metric,
                     nmi = .nmi,
                     ari = .ari
      )
      clustering.scores <- sapply(clustering.results, FUN = FUN, cell.var.vec,
                                  simplify = "numeric", USE.NAMES = TRUE)
      best.idx <- which.max(clustering.scores)
      best.score <- clustering.scores[best.idx]
      res <- names(best.score)
      best.res <- sub("^res.", "", res)

      message(sprintf("\tResolution %s maximized %s (%.2f)\n", best.res,
                      toupper(metric), best.score)[verbose],
              appendLF = FALSE)

      clustering.results %>%
        select({{ best.idx }}) %>% rename(!!englue(cluster.name) := colnames(.))
      # names(best.clustering) <- paste0(names(best.clustering), "_", metric)
      # best.clustering
      ### end
    }, simplify = FALSE, USE.NAMES = TRUE) %>% bind_cols()
  }) %>% bind_cols()

  # FUN <- switch (optimisation.metric,
  #                nmi = SeuratIntegrate:::.nmi,
  #                ari = SeuratIntegrate:::.ari
  # )
  # clustering.scores <- sapply(clustering.results, FUN = FUN, cell.var.vec,
  #                             simplify = "numeric", USE.NAMES = TRUE)
  # best.idx <- which.max(clustering.scores)
  # best.score <- clustering.scores[best.idx]
  # best.res <- sub("^res.", "", names(best.score))
  #
  # message(sprintf("Resolution %s maximized %s (%.2f)\n", best.res,
  #                 toupper(optimisation.metric), best.score)[verbose],
  #         appendLF = FALSE)
  #
  # clustering.results <- clustering.results[, best.idx, drop = FALSE]
  # ### end
  #
  # cluster.name <- cluster.name %||% paste(graph.name, names(x = clustering.results),
  #                                         sep = "_")
  # names(x = clustering.results) <- cluster.name
  # idents.use <- names(x = clustering.results)[ncol(x = clustering.results)]
  # object[[]] <- clustering.results
  # Idents(object = object, replace = TRUE) <- object[[idents.use,
  #                                                    drop = TRUE]]
  # levels <- levels(x = object)
  # levels <- tryCatch(expr = as.numeric(x = levels), warning = function(...) {
  #   return(levels)
  # }, error = function(...) {
  #   return(levels)
  # })
  # Idents(object = object) <- factor(x = Idents(object = object),
  #                                   levels = sort(x = levels))
  # object[["seurat_clusters"]] <- Idents(object = object)
  object <- AddMetaData(object = object, metadata = best.clusters)
  cmd <- LogSeuratCommand(object = object, return.command = TRUE)
  slot(object = cmd, name = "assay.used") <- DefaultAssay(object = object[[graph.name]])
  object[[slot(object = cmd, name = "name")]] <- cmd
  return(object)
}
