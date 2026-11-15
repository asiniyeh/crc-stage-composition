# =============================================================================
# 00_config.R  -  settings shared by every script
# Manuscript 1: stage-associated gene expression in colorectal cancer,
# adjusted for tumour purity and consensus molecular subtype (CMS)
# =============================================================================
# Every threshold below was fixed in the study protocol BEFORE looking at data.
# If you change one, record the change and the reason in the manuscript.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

# Bioconductor packages loaded later (matrixStats, AnnotationDbi, S4Vectors) mask some dplyr verbs.
# Bind the dplyr versions explicitly so every script gets the right one.
count <- dplyr::count; select <- dplyr::select; filter <- dplyr::filter
rename <- dplyr::rename; slice <- dplyr::slice; desc <- dplyr::desc

# ---- folders ---------------------------------------------------------------
PROJECT_DIR <- normalizePath(".", mustWork = TRUE)   # run scripts from the project folder
DIR_INPUT   <- file.path(PROJECT_DIR, "input")
DIR_RAW     <- file.path(PROJECT_DIR, "data_raw")    # downloads (large, not shared)
DIR_DATA    <- file.path(PROJECT_DIR, "data")        # cleaned R objects
DIR_RESULTS <- file.path(PROJECT_DIR, "results")     # tables to send back
DIR_FIGS    <- file.path(PROJECT_DIR, "figures")     # figures to send back
for (d in c(DIR_RAW, DIR_DATA, DIR_RESULTS, DIR_FIGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- cohorts ---------------------------------------------------------------
TCGA_PROJECTS <- c("TCGA-COAD", "TCGA-READ")
GEO_COHORTS   <- c("GSE39582", "GSE17536")

# ---- pre-specified thresholds ----------------------------------------------
FDR_CUTOFF        <- 0.05   # adjusted p-value for differential expression
LFC_CUTOFF        <- 0.5    # minimum |log2 fold change|, late vs early stage
STROMA_RHO_MAX    <- 0.4    # panel genes must have |Spearman rho| with stromal score below this
MIN_COUNT         <- 10     # gene kept if it has >= MIN_COUNT reads ...
MIN_SAMPLE_FRAC   <- 0.10   # ... in at least this fraction of tumours
VALIDATION_P      <- 0.05   # nominal p in a validation cohort, same direction = replicated
GSEA_MIN_SIZE     <- 15
GSEA_MAX_SIZE     <- 500
CV_FOLDS          <- 5      # outer and inner folds for the LASSO classifier
SEED              <- 20261114

set.seed(SEED)

# ---- helpers used by several scripts ---------------------------------------

# "Stage IIIB" -> "III";  "stage 2" / "2" -> "II";  anything else -> NA
parse_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^STAGE\\s*", "", x)
  x <- sub("^AJCC\\s*", "", x)
  out <- rep(NA_character_, length(x))
  num <- suppressWarnings(as.integer(sub("[^0-9].*$", "", x)))
  roman <- sub("[ABC0-9]*$", "", x)
  out[roman %in% c("I", "II", "III", "IV")] <- roman[roman %in% c("I", "II", "III", "IV")]
  out[is.na(out) & num %in% 1:4] <- c("I", "II", "III", "IV")[num[is.na(out) & num %in% 1:4]]
  out
}

stage_group <- function(stage4) {
  factor(ifelse(stage4 %in% c("I", "II"), "early",
         ifelse(stage4 %in% c("III", "IV"), "late", NA)),
         levels = c("early", "late"))
}

# Turn event columns written as text ("death", "no death", "1", "0", "relapse") into 0/1
to_event <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- suppressWarnings(as.numeric(x))
  txt <- is.na(out) & !is.na(x) & x != "" & x != "na" & x != "n/a"
  out[txt] <- ifelse(grepl("^no|^none|censor|alive|living|^0", x[txt]), 0, 1)
  out
}

# Find one phenotype column by regular expression; stop with a helpful list if absent
find_col <- function(df, pattern, required = TRUE) {
  hit <- grep(pattern, colnames(df), ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0) {
    if (required) {
      stop("No column matches '", pattern, "'. Available columns:\n  ",
           paste(colnames(df), collapse = "\n  "),
           "\nEdit the pattern in the script and run it again.", call. = FALSE)
    }
    return(NA_character_)
  }
  if (length(hit) > 1) message("  pattern '", pattern, "' matched several columns; using '", hit[1], "'")
  hit[1]
}

# Collapse a probe/feature matrix to one row per gene symbol: keep the row with the highest mean
collapse_by_symbol <- function(mat, symbols) {
  keep <- !is.na(symbols) & symbols != "" & !grepl("///", symbols)
  mat <- mat[keep, , drop = FALSE]; symbols <- symbols[keep]
  ord <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
  mat <- mat[ord, , drop = FALSE]; symbols <- symbols[ord]
  first <- !duplicated(symbols)
  out <- mat[first, , drop = FALSE]
  rownames(out) <- symbols[first]
  out
}

save_table <- function(df, name) {
  path <- file.path(DIR_RESULTS, paste0(name, ".csv"))
  readr::write_csv(df, path)
  message("  saved ", path)
  invisible(path)
}

save_fig <- function(plot, name, width = 7, height = 5) {
  path <- file.path(DIR_FIGS, paste0(name, ".pdf"))
  ggplot2::ggsave(path, plot, width = width, height = height)
  ggplot2::ggsave(sub("\\.pdf$", ".png", path), plot, width = width, height = height, dpi = 300)
  message("  saved ", path, " (+ .png)")
  invisible(path)
}

theme_set(theme_bw(base_size = 11))
