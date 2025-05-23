% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conda_S4_helpers.R
\name{resetCache}
\alias{resetCache}
\title{Reset a conda environment}
\usage{
resetCache(method = known.methods)
}
\arguments{
\item{method}{the name of the method to update. One of "combat", "harmony",
"mnn", "bbknn", "scvi", "scanvi", "scanorama", "trvae"}
}
\value{
a \linkS4class{CondaEnvManager} object (invisibly)
}
\description{
Reset (unset) a conda environment linked to a provided method. Update cache
on disk and in package environment
}
