# =============================================================================
# 01_download_clean.R  -  download TCGA + GEO data and build one clean table per cohort
# Output: data/tcga.rds, data/GSE39582.rds, data/GSE17536.rds,
#         results/table1_characteristics.csv, results/flow_counts.csv
# Time: 30-90 minutes (TCGA download is about 2-3 GB)
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(GEOquery)
  library(Biobase)
})

flow <- list()   # sample counts at each step, for the study-flow figure
log_flow <- function(cohort, step, n) {
  flow[[length(flow) + 1]] <<- tibble(cohort = cohort, step = step, n = n)
  message(sprintf("  [%s] %-45s n = %d", cohort, step, n))
}

# =============================================================================
# 1. TCGA-COAD + TCGA-READ  (RNA-seq STAR counts)
# =============================================================================
# GDCprepare() loads all six expression measures for every file at once and needs
# >8 GB RAM for ~700 files. Instead we read only the two columns we use
# (raw counts, TPM) from each downloaded file, which needs <2 GB.
message("\n== TCGA ==")
read_star <- function(path) {
  d <- data.table::fread(path, skip = 1, select = c("gene_id", "gene_name", "gene_type", "unstranded", "tpm_unstranded"),
                         showProgress = FALSE)
  d[grepl("^ENSG", d$gene_id), ]
}
files <- list(); clin_raw <- list()
for (proj in TCGA_PROJECTS) {
  q <- GDCquery(
    project       = proj,
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  GDCdownload(q, directory = DIR_RAW, files.per.chunk = 10)   # skips files already present
  r <- getResults(q)
  r$path <- file.path(DIR_RAW, proj, "Transcriptome_Profiling", "Gene_Expression_Quantification", r$id, r$file_name)
  r$project <- proj
  files[[proj]] <- r[, c("id", "cases", "project", "path")]
  clin_raw[[proj]] <- GDCquery_clinic(proj, type = "clinical")
}
files <- do.call(rbind, files)
stopifnot(all(file.exists(files$path)))
message("  reading ", nrow(files), " STAR count files")
first <- read_star(files$path[1])
genes <- data.frame(gene_id = first$gene_id, gene_name = first$gene_name, gene_type = first$gene_type)
genes$ensg <- sub("\\..*$", "", genes$gene_id)
counts <- matrix(0L, nrow(genes), nrow(files), dimnames = list(genes$ensg, files$cases))
tpm    <- matrix(0,  nrow(genes), nrow(files), dimnames = list(genes$ensg, files$cases))
for (i in seq_len(nrow(files))) {
  d <- read_star(files$path[i])
  stopifnot(identical(d$gene_id, genes$gene_id))
  counts[, i] <- as.integer(d$unstranded)
  tpm[, i]    <- d$tpm_unstranded
  if (i %% 100 == 0) message("    ", i, " files read")
}

# clinical data are per patient
cl <- dplyr::bind_rows(lapply(clin_raw, function(x) x[, intersect(c("submitter_id", "ajcc_pathologic_stage", "vital_status",
  "days_to_death", "days_to_last_follow_up", "age_at_index", "sex_at_birth", "gender", "tissue_or_organ_of_origin"), colnames(x))]))
cl <- cl[!duplicated(cl$submitter_id), ]
patient <- substr(files$cases, 1, 12)
cd <- cl[match(patient, cl$submitter_id), ]
sex_col <- intersect(c("sex_at_birth", "gender"), colnames(cd))[1]

# MSI status: published TCGA subtype table (TCGA Network 2012); covers only part of the cohort
st <- tryCatch(TCGAquery_subtype("coad"), error = function(e) NULL)
msi_raw <- if (!is.null(st)) as.character(st$MSI_status[match(patient, substr(st$patient, 1, 12))]) else rep(NA_character_, length(patient))
message("  TCGA samples with MSI status: ", sum(!is.na(msi_raw) & msi_raw %in% c("MSI-H", "MSI-L", "MSS")), " of ", length(patient))

side_map <- function(site) {
  s <- tolower(site)
  ifelse(grepl("cecum|ascending|hepatic|transverse", s), "right",
  ifelse(grepl("splenic|descending|sigmoid|rect", s), "left", NA))
}

clin_all <- tibble(
  sample_id  = files$cases,
  patient    = patient,
  sample_code = substr(files$cases, 14, 15),
  project    = files$project,
  stage4     = parse_stage(cd$ajcc_pathologic_stage),
  age        = as.numeric(cd$age_at_index),
  sex        = tolower(cd[[sex_col]]),
  organ      = cd$tissue_or_organ_of_origin,
  side       = side_map(cd$tissue_or_organ_of_origin),
  site       = ifelse(grepl("rect", tolower(cd$tissue_or_organ_of_origin)) | files$project == "TCGA-READ", "rectum", "colon"),
  msi        = case_when(msi_raw == "MSI-H" ~ "MSI-H",
                         msi_raw %in% c("MSS", "MSI-L") ~ "MSS",
                         TRUE ~ NA_character_),
  os_time    = as.numeric(ifelse(is.na(cd$days_to_death), cd$days_to_last_follow_up, cd$days_to_death)) / 30.44,
  os_event   = as.numeric(tolower(cd$vital_status) == "dead")
)
log_flow("TCGA", "All RNA-seq samples downloaded", nrow(clin_all))

# normal mucosa (sample code 11) is kept aside for script 09 (tumour vs normal)
is_normal <- clin_all$sample_code == "11"
normal <- list(counts = counts[, is_normal, drop = FALSE], clin = clin_all[is_normal, ])

# primary tumours only, one per patient, stage known
tum_idx <- which(clin_all$sample_code == "01")
log_flow("TCGA", "Primary tumour samples", length(tum_idx))
tum_idx <- tum_idx[order(clin_all$sample_id[tum_idx])]          # "01A" before "01B"
tum_idx <- tum_idx[!duplicated(clin_all$patient[tum_idx])]      # first aliquot per patient
keep <- seq_len(nrow(clin_all)) %in% tum_idx
log_flow("TCGA", "One sample per patient", sum(keep))
keep <- keep & !is.na(clin_all$stage4)
log_flow("TCGA", "Stage I-IV recorded", sum(keep))

clin <- clin_all[keep, ] %>% mutate(cohort = "TCGA", stage_group = stage_group(stage4))
tcga <- list(
  counts = counts[, keep], tpm = tpm[, keep], genes = genes, clin = clin, normal = normal
)

# published ABSOLUTE purity (Aran et al. 2015), used only in a sensitivity analysis
tcga$purity_published <- tryCatch({
  data("Tumor.purity", package = "TCGAbiolinks", envir = environment())
  tp <- Tumor.purity
  tp$sample16 <- substr(tp$Sample.ID, 1, 16)
  tibble(sample_id = clin$sample_id,
         absolute = suppressWarnings(as.numeric(sub(",", ".", tp$ABSOLUTE[match(substr(clin$sample_id, 1, 16), tp$sample16)]))))
}, error = function(e) { message("  published purity table not available: ", conditionMessage(e)); NULL })

saveRDS(tcga, file.path(DIR_DATA, "tcga.rds"))
rm(counts, tpm); gc()

# =============================================================================
# 2. GEO cohorts (Affymetrix U133 Plus 2.0, RMA-normalised series matrix)
# =============================================================================
# Files are fetched from NCBI's download server (ftp.ncbi.nlm.nih.gov) so the analysis
# does not depend on the GEO web interface, and probe -> gene symbols come from the
# official GPL570 annotation file.
geo_dir <- file.path(DIR_RAW, "geo"); dir.create(geo_dir, showWarnings = FALSE)
geo_get <- function(rel_path) {
  dest <- file.path(geo_dir, basename(rel_path))
  if (!file.exists(dest)) download.file(paste0("https://ftp.ncbi.nlm.nih.gov/geo/", rel_path), dest, mode = "wb")
  dest
}
gpl_file <- geo_get("platforms/GPLnnn/GPL570/annot/GPL570.annot.gz")
gpl_lines <- readLines(gzfile(gpl_file), n = 200)
gpl <- read.delim(gzfile(gpl_file), skip = grep("^ID\t", gpl_lines)[1] - 1, quote = "",
                  comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
gpl <- gpl[!grepl("^!", gpl$ID), ]
probe2symbol <- setNames(gpl[["Gene symbol"]], gpl$ID)

load_geo <- function(gse) {
  message("\n== ", gse, " ==")
  stub <- paste0(substr(gse, 1, nchar(gse) - 3), "nnn")
  f <- geo_get(sprintf("series/%s/%s/matrix/%s_series_matrix.txt.gz", stub, gse, gse))
  eset <- getGEO(filename = f, getGPL = FALSE)
  ex <- exprs(eset)
  if (max(ex, na.rm = TRUE) > 50) ex <- log2(ex + 1)       # make sure values are log2
  expr <- collapse_by_symbol(ex, unname(probe2symbol[rownames(ex)]))
  message("  genes after collapsing probes: ", nrow(expr))
  pd <- pData(eset)
  message("  phenotype columns available:\n    ", paste(grep(":ch1$", colnames(pd), value = TRUE), collapse = "\n    "))
  list(expr = expr, pd = pd)
}

# ---- GSE39582 (Marisa et al. 2013) -----------------------------------------
g <- load_geo("GSE39582")
pd <- g$pd
log_flow("GSE39582", "All samples", nrow(pd))
# the 19 non-tumoral mucosa samples carry their patient's stage, so filter on the sample source
is_tumour <- !grepl("non tumoral", tolower(pd[[find_col(pd, "^source_name_ch1$")]]))
loc <- tolower(pd[[find_col(pd, "tumor\\.location:ch1")]])
mmr <- tolower(pd[[find_col(pd, "mmr\\.status:ch1")]])
chemo_col <- find_col(pd, "chemotherapy\\.adjuvant:ch1", required = FALSE)
clin <- tibble(
  sample_id   = rownames(pd),
  cohort      = "GSE39582",
  stage4      = parse_stage(pd[[find_col(pd, "tnm\\.stage:ch1")]]),
  age         = suppressWarnings(as.numeric(pd[[find_col(pd, "age\\.at\\.diagnosis")]])),
  sex         = tolower(pd[[find_col(pd, "^sex:ch1")]]),
  side        = ifelse(loc == "proximal", "right", ifelse(loc == "distal", "left", NA)),
  site        = "colon",
  msi         = ifelse(grepl("dmmr", mmr), "MSI-H", ifelse(grepl("pmmr", mmr), "MSS", NA)),
  chemo       = if (is.na(chemo_col)) NA_character_ else
                  dplyr::recode(tolower(pd[[chemo_col]]), y = "yes", n = "no", .default = NA_character_),
  os_time     = suppressWarnings(as.numeric(pd[[find_col(pd, "os\\.delay")]])),
  os_event    = to_event(pd[[find_col(pd, "os\\.event")]]),
  rfs_time    = suppressWarnings(as.numeric(pd[[find_col(pd, "rfs\\.delay")]])),
  rfs_event   = to_event(pd[[find_col(pd, "rfs\\.event")]])
)
keep <- is_tumour
log_flow("GSE39582", "Tumour samples", sum(keep))
keep <- keep & !is.na(clin$stage4)
log_flow("GSE39582", "Stage I-IV recorded", sum(keep))
clin <- clin[keep, ] %>% mutate(stage_group = stage_group(stage4))
saveRDS(list(expr = g$expr[, clin$sample_id], clin = clin), file.path(DIR_DATA, "GSE39582.rds"))

# ---- GSE17536 (Smith et al. 2010, Moffitt) ---------------------------------
g <- load_geo("GSE17536")
pd <- g$pd
log_flow("GSE17536", "All samples", nrow(pd))
clin <- tibble(
  sample_id = rownames(pd),
  cohort    = "GSE17536",
  stage4    = parse_stage(pd[[find_col(pd, "ajcc")]]),
  age       = suppressWarnings(as.numeric(pd[[find_col(pd, "^age")]])),
  sex       = tolower(pd[[find_col(pd, "gender|^sex")]]),
  side      = NA_character_,
  site      = NA_character_,
  msi       = NA_character_,
  chemo     = NA_character_,
  os_time   = suppressWarnings(as.numeric(pd[[find_col(pd, "overall survival follow-up time")]])),
  os_event  = to_event(pd[[find_col(pd, "overall_event")]]),
  rfs_time  = suppressWarnings(as.numeric(pd[[find_col(pd, "dfs_time")]])),
  rfs_event = to_event(pd[[find_col(pd, "dfs_event")]])
)
keep <- !is.na(clin$stage4)
log_flow("GSE17536", "Stage I-IV recorded", sum(keep))
clin <- clin[keep, ] %>% mutate(stage_group = stage_group(stage4))
saveRDS(list(expr = g$expr[, clin$sample_id], clin = clin), file.path(DIR_DATA, "GSE17536.rds"))

# =============================================================================
# 3. Table 1 and flow counts
# =============================================================================
all_clin <- bind_rows(
  readRDS(file.path(DIR_DATA, "tcga.rds"))$clin %>% mutate(chemo = NA_character_, rfs_time = NA, rfs_event = NA),
  readRDS(file.path(DIR_DATA, "GSE39582.rds"))$clin,
  readRDS(file.path(DIR_DATA, "GSE17536.rds"))$clin
)
pct <- function(x, level) sprintf("%d (%.0f%%)", sum(x == level, na.rm = TRUE), 100 * mean(x == level, na.rm = TRUE))
table1 <- all_clin %>%
  group_by(cohort) %>%
  summarise(
    n                 = n(),
    age_median_IQR    = sprintf("%.0f (%.0f-%.0f)", median(age, na.rm = TRUE), quantile(age, .25, na.rm = TRUE), quantile(age, .75, na.rm = TRUE)),
    female            = pct(sex, "female"),
    stage_I           = pct(stage4, "I"),
    stage_II          = pct(stage4, "II"),
    stage_III         = pct(stage4, "III"),
    stage_IV          = pct(stage4, "IV"),
    MSI_H             = pct(msi, "MSI-H"),
    MSI_missing       = sum(is.na(msi)),
    right_sided       = pct(side, "right"),
    rectum            = pct(site, "rectum"),
    os_events         = sum(os_event, na.rm = TRUE),
    rfs_events        = sum(rfs_event, na.rm = TRUE),
    followup_median_m = round(median(os_time, na.rm = TRUE), 1)
  )
save_table(table1, "table1_characteristics")
save_table(bind_rows(flow), "flow_counts")

# quick check the sender can read: event coding per cohort
print(all_clin %>% dplyr::count(cohort, os_event))
print(all_clin %>% dplyr::count(cohort, rfs_event))
message("\n01 done. Check the event counts above look sensible (0 = censored, 1 = event).")
