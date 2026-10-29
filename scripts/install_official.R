# Optional local R dependencies for the original-code reproduction.
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(file_arg)) stop("Run with Rscript scripts/install_official.R")
script <- normalizePath(sub("^--file=", "", file_arg[[1]]))
root <- dirname(dirname(script))
lib <- file.path(root, ".cache", "r-library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
cran <- c("BiocManager", "data.table", "dplyr", "tidyr", "tibble", "R.utils",
          "ROCR", "caTools", "yardstick", "yaml")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, lib = lib, type = .Platform$pkgType)
bioc <- c("GenomicRanges", "BiocParallel")
missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) BiocManager::install(missing, lib = lib, ask = FALSE,
                                       update = FALSE, type = .Platform$pkgType)
cat("Local library:", lib, "\nNext: Rscript scripts/reproduce_official.R\n")
