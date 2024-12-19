# Seurat-Integrate
R package gathering a set of wrappers to apply various integration methods to Seurat objects and rate integration obtained with such methods

## Installation
### 1) Remotely from github with R
Install SeuratIntegrate from github directly:
```R
install.packages(c("remotes", "BiocManager"))
remotes::install_github("cbib/Seurat-Integrate", dependancies = NA, repos = BiocManager::repositories()) 
```
If you prefer, or if it does not work, you can use the alternative way described below.
### 2) From a local copy using the remotes R package
The installation process encompasses two steps, namely:
 - Clone or download the repository or the [latest release](https://github.com/cbib/Seurat-Integrate/releases/tag/0.3.1)
 - Install the package

First off, go to `some/path/to/download/the/repository`. Then, clone or download the repository or a release.
```shell
# Example: clone
git clone git@github.com:cbib/Seurat-Integrate.git
```
Depending on what you chose to download, you can either obtain the source code in a `Seurat-Integrate/` folder or a compressed version ending with `.targ.gz`. Regardless, you can use the command below:
```R
install.packages(c("remotes", "BiocManager"))
path <- 'Seurat-Integrate/'   # or SeuratIntegrate_X.X.X.tar.gz (X.X.X being the version)
remotes::install_local(path = path, dependencies = NA, repos = BiocManager::repositories())
```
Done ! Open a R session and try it out !
