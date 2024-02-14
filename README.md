# Seurat-Integrate
R package gathering a set of wrappers to apply various integration methods to Seurat objects (and rate such methods)

## Installation
The installation process encompasses three steps, namely:
 - Clone or download the repository
 - Build the package using R
 - Install the package

First off, go to `some/path/to/download/the/repository`. Then,
```shell
# Clone
git clone git@github.com:fspecque/Seurat-Integrate.git

# Build and install (assuming that the R binary is in your PATH)
R CMD build Seurat-Integrate
version=0.1.0  # current version
R CMD INSTALL SeuratIntegrate_${version}.tar.gz
```
Done ! Open a R session and try it out !

## Development
To particiape in developing this package, you need to get farmiliar with a few principles. Don't hesitate to have a look at [this book](https://r-pkgs.org/), it is a good information resource.

> You can check and modify the [TODO issue](/../../issues/1) 

### New functions
For now, any new function must be placed inside its own `.R` file inside the `R/` folder.

To document new functions, [roxygen2](https://roxygen2.r-lib.org) provide useful `@tags` that are place above the function definition to ease the generation of the manuals (`?myfunction`). `R/SeuratIntegrate-package.R` provides default parameter values shared between integration methods and required by `Seurat::IntegrateLayers`.

**Be sure to place `#' @inheritParams integration-method` in the roxygen2 documentation block** before the function definition in the R file to include them.

### Update files
During the development process, when you working directory is `Seurat-Integrate/`, **be sure to run those lines** inside a R session.
```shell
roxygen2::roxygenise()
devtools::document()
```
> Those two functions update the `NAMESPACE` file and the functions manual in `man/`

When you push your new function (after running the last two commands), be sure to push what changed in `R/`, `man/`, `NAMESPACE` (and `DESCRIPTION` if applicable)
