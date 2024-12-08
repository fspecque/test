#' @importFrom utils globalVariables

if (getRversion() >= "2.15.1")
  globalVariables(c("n_cells", "n_batches", "k0_bound", ".",
                    "df.mtdt", "idcol", "dimred", "dimvar",
                    "method", "var", "H", "proba", "sigmas", "rhos",
                    "rows", "cols", "vals", "distances",
                    "Integration", "Overall.score", "Score", "y", "Score.type",
                    "RECT_x", "RECT_y", "sub.object", "df.score"))
