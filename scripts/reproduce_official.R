# Execute original CRISPR_comparison example numerical routines, unchanged.
# Snakemake orchestration is replaced by this small context; upstream R code is not.
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(file_arg)) stop("Run with Rscript scripts/reproduce_official.R")
script <- normalizePath(sub("^--file=", "", file_arg[[1]]))
root <- dirname(dirname(script))
cache <- file.path(root, ".cache")
lib <- file.path(cache, "r-library")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
required <- c("data.table", "dplyr", "tidyr", "tibble", "R.utils", "ROCR",
              "caTools", "yardstick", "yaml", "GenomicRanges", "BiocParallel")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R dependencies: ", paste(missing, collapse = ", "),
       "\nInstall locally with: Rscript scripts/install_official.R", call. = FALSE)
}
upstream <- file.path(cache, "CRISPR_comparison")
commit <- "50587422e6b11259ead6fbc6f867681c788f39b7"
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
git <- Sys.which("git")
if (!nzchar(git)) stop("git must be installed and on PATH.")
if (!dir.exists(file.path(upstream, ".git"))) {
  if (dir.exists(upstream)) stop("Cache exists but is not a Git checkout: ", upstream)
  status <- system2(git, c("clone", "--quiet", "https://github.com/EngreitzLab/CRISPR_comparison.git", shQuote(upstream)))
  if (status != 0) stop("Could not clone the public upstream repository.")
  status <- system2(git, c("-C", shQuote(upstream), "checkout", "--quiet", commit))
  if (status != 0) stop("Could not check out the pinned upstream commit.")
}
observed_commit <- system2(git, c("-C", shQuote(upstream), "rev-parse", "HEAD"), stdout = TRUE)
if (!identical(observed_commit, commit)) {
  stop("Existing upstream cache is not the pinned commit. Inspect it before changing it.")
}
changed <- system2(git, c("-C", shQuote(upstream), "status", "--porcelain", "--untracked-files=no"), stdout = TRUE)
if (length(changed)) stop("Tracked upstream files are modified; refusing to label this an exact reproduction.")
result_dir <- file.path(root, "results", "official")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TZ = "UTC")
set.seed(20261029)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ROCR)
  library(caTools)
  library(GenomicRanges)
})
setwd(upstream)
config <- yaml::read_yaml("config/config.yml")
ex <- config$comparisons$example
methods::setClass("OfficialExampleContext", slots = c(
  input = "list", output = "list", params = "list", log = "list",
  config = "list", wildcards = "list", scriptdir = "character"))
snakemake <- methods::new("OfficialExampleContext",
  input = list(predictions = unname(ex$pred), experiment = ex$expt,
    tss_universe = ex$tss_universe, gene_universe = ex$gene_universe,
    pred_config = ex$pred_config, cell_type_mapping = unname(ex$cell_type_mapping),
    expressed_genes = ex$expressed_genes),
  output = list(merged = file.path(result_dir, "expt_pred_merged.txt.gz")),
  params = list(pos_col = "Regulated", include_col = NULL,
    pred_format = "ENCODE", filter_tss = TRUE, filter_include_col = FALSE),
  log = list(file.path(result_dir, "mergePredictionsWithExperiment.log")),
  config = config, wildcards = list(comparison = "example"),
  scriptdir = file.path(upstream, "workflow/scripts"))
source(file.path(snakemake@scriptdir, "mergePredictionsWithExperiment.R"))

# The example has no optional annotations, so that upstream rule merely copies
# its input. Read it directly. Load original plotting helper definitions without
# unused plot-package imports: no function bodies are modified.
plot_expr <- parse(file.path(snakemake@scriptdir, "crisprComparisonPlotFunctions.R"))
for (expr in plot_expr) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
    eval(expr, envir = .GlobalEnv)
  }
}
source(file.path(snakemake@scriptdir, "crisprComparisonBootstrapFunctions.R"))
merged <- fread(snakemake@output$merged, colClasses = c(ValidConnection = "character"))
pred_config <- importPredConfig(ex$pred_config,
  expr = "baseline.nearestExprTSS" %in% merged$pred_uid, include_col = "include")
processed <- processMergedData(merged, pred_config,
  filter_valid_connections = TRUE, include_missing_predictions = TRUE)
pred_config <- subset(pred_config, pred_uid %in% unique(processed$pred_uid))
fwrite(processed, file.path(result_dir, "processed_merged.tsv.gz"), sep = "\t", na = "NA")
pr <- applyCellTypes(processed, .fun = calcPRCurves, pred_config = pred_config, pos_col = "Regulated")
fwrite(rbindlist(pr, idcol = "ExperimentCellType"), file.path(result_dir, "pr_table.tsv"), sep = "\t", na = "NA")
perf <- applyCellTypes(processed, .fun = makePRSummaryTableBS,
  pred_config = pred_config, pos_col = "Regulated", R = 1000)
stopifnot(length(perf) > 0)
perf <- rbindlist(perf, idcol = "cell_type")
fwrite(perf, file.path(result_dir, "performance_summary.tsv"), sep = "\t", na = "NA")
counts <- processed[, .(n = .N, positives = sum(Regulated),
  negatives = sum(!Regulated), missing_predictions = sum(Prediction == 0)),
  by = .(ExperimentCellType, pred_uid)]
fwrite(counts, file.path(result_dir, "counts.tsv"), sep = "\t")
stopifnot(all(counts$n == 3852L), all(counts$positives == 101L))
stopifnot(abs(perf[pred_uid == "ABC.powerlaw.Score", AUPRC] - 0.588235447008955) < 1e-10)
stopifnot(abs(perf[pred_uid == "baseline.distToTSS", AUPRC] - 0.401794369883952) < 1e-10)
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
loaded <- unique(c(required, loadedNamespaces()))
versions <- data.frame(Package = loaded,
  Version = vapply(loaded, function(x) as.character(packageVersion(x)), character(1)))
fwrite(versions, file.path(result_dir, "loaded_versions.tsv"), sep = "\t")
writeLines(c(paste("upstream_commit", commit, sep = "\t"),
  "seed\t20261029", "bootstrap_replicates\t1000",
  "method\toriginal merge script and original numerical report functions; no source edits",
  "html_report_rendered\tFALSE"), file.path(result_dir, "provenance.tsv"))
print(counts)
print(perf)
cat("Official numerical reproduction verified. Outputs: results/official/\n")
