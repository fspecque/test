################################################################################
##########################                            ##########################
#######                          useful objects                          #######
###                                                                          ###
#                                                                              #

known.methods <- c("combat", "harmony", "mnn", "bbknn", "scvi", "scanvi",
                   "scanorama", "trvae")
needs.conda <- sapply(known.methods, function(method) ifelse(
  method %in% c("scanorama", "bbknn", "scvi", "scanvi", "trvae"), T, F
))

.conda_requirements <- list(
  bbknn = list(
    packages = c("python", "scipy", "scanpy", "bbknn"),
    channels = c("conda-forge", "bioconda"),
    pip = c()
  ),
  scvi = list(
    packages = c("python", "scipy", "scanpy", "scvi-tools", "mkl",
                 "mkl-include", "setuptools", "cmake"),
    channels = c("conda-forge"),
    pip = c()
  ),
  scanorama = list(
    packages = c("python", "scanorama", "scipy"),
    channels = c("conda-forge", "bioconda"),
    pip = c()
  ),
  trvae = list(
    packages = c("python", "pip", "numpy", "pytorch", "torchaudio",
                 "torchvision", "pytorch-cuda", "scvi-tools"),
    channels = c("conda-forge", "bioconda", "pytorch", "nvidia"),
    pip = c("scarches")
  )
)
.conda_requirements$scanvi <- .conda_requirements$scvi

.conda_default_envname <- list(
  bbknn = "SeuratIntegrate_bbknn",
  scvi = "SeuratIntegrate_scvi-tools",
  scanorama = "SeuratIntegrate_scanorama",
  trvae = "SeuratIntegrate_trvae"
)
.conda_default_envname$scanvi <- .conda_default_envname$scvi
################################################################################
##########################                            ##########################
#######                         class definitions                        #######
###                                                                          ###
#                                                                              #

setClass(
  Class = "CondaEnvSlot",
  slots = list(
    value = 'character',
    valid = 'logical'
  )
)

#' Encapsulates information about a conda environment
#'
#' @description
#' The \code{CondaEnv} class provides a basic structure to store information
#' about a conda environment for a given method. It is designed to be used by
#' the \link[=CondaEnvManager-class]{CondaEnvManager} class.
#'
#' @slot method 1-length character. Indicating the method this \code{CondaEnv}
#' refers to
#' @slot needs.conda logical. Whether the \code{method} needs conda (yes if
#' python-based, no if R-based)
#' @slot conda.bin single string. path to conda binary. If empty (\code{""}),
#' \code{NULL} or \code{"conda"}, \code{"auto"} is passed to
#' \link[reticulate]{conda_binary} to find path to conda binary in \code{PATH}
#' @slot conda.env.name 1-length character. Name of the conda environment
#' @slot conda.env.path 1-length character. Path to the conda environment
#' @slot is.valid logical. Whether the environment is valid. Should not be set
#' by the user. Set by \link{checkCondaEnv}.
#' @exportClass CondaEnv
#' @keywords internal
#' @seealso \link[=CondaEnv]{CondaEnv()}, \linkS4class{CondaEnvManager},
#' \link[reticulate:conda-tools]{conda_binary()}
setClass(
  Class = "CondaEnv",
  slots = list(
    method = 'character',
    needs.conda = 'logical',
    conda.bin = 'CondaEnvSlot',
    conda.env.name = 'CondaEnvSlot',
    conda.env.path = 'CondaEnvSlot',
    is.valid = 'logical'
  )
)

#' Manager of conda environments for python-based integration methods
#'
#' @description
#' The \code{CondaEnvManager} class provides a handy way to set up, store and use
#' conda environments for python-based integration methods. It is designed to be
#' set up and modified via helper functions, not directly by the user.
#'
#' @slot combat CondaEnv. For R-based combat method, nothing to set up.
#' @slot harmony CondaEnv. For R-based harmony method, nothing to set up.
#' @slot mnn CondaEnv. For R-based MNN method, nothing to set up.
#' @slot bbknn CondaEnv. For python-based bbknn method.
#' @slot scvi CondaEnv. For python-based SCVI method. Can be shared with SCANVI.
#' @slot scanvi CondaEnv. For python-based SCANVI method. Can be shared with SCVI.
#' @slot scanorama CondaEnv. For python-based Scanorama method
#' @slot trvae CondaEnv. For python-based trVAE method
#' @exportClass CondaEnvManager
#' @keywords internal
#' @seealso \linkS4class{CondaEnv}
setClass(
  Class = "CondaEnvManager",
  slots = sapply(known.methods, function(whatever) 'CondaEnv', simplify = F)
)
