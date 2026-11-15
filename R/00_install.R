# =============================================================================
# 00_install.R  -  run ONCE to install every package the pipeline needs
# Needs R >= 4.3 and an internet connection. Takes 20-40 minutes the first time.
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "dplyr", "tibble", "readr", "tidyr", "stringr", "ggplot2", "ggrepel", "patchwork",
  "survival", "survminer", "glmnet", "pROC", "metafor", "msigdbr", "tidyestimate",
  "remotes", "renv"
)
bioc_pkgs <- c(
  "TCGAbiolinks", "SummarizedExperiment", "GEOquery", "Biobase", "DESeq2", "apeglm",
  "limma", "edgeR", "fgsea", "BiocParallel"
)

new_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(new_cran)) install.packages(new_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
new_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(new_bioc)) BiocManager::install(new_bioc, update = FALSE, ask = FALSE)

# CMScaller is on GitHub only (Eide et al. 2017, Sci Rep)
if (!requireNamespace("CMScaller", quietly = TRUE)) {
  remotes::install_github("Lothelab/CMScaller", upgrade = "never")
}

# Check everything loads
all_pkgs <- c(cran_pkgs, bioc_pkgs, "CMScaller")
ok <- vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (all(ok)) {
  message("\nAll ", length(ok), " packages installed.")
} else {
  message("\nThese packages failed to install: ", paste(names(ok)[!ok], collapse = ", "))
  message("Copy the error messages above and send them to Claude.")
}

# Record exact versions for the manuscript (reproducibility)
writeLines(capture.output(sessionInfo()), "session_info_install.txt")
