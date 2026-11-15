# =============================================================================
# 02_purity_cms.R  -  stromal/immune scores (ESTIMATE) and consensus molecular subtype
# Output: adds stromal, immune, estimate, purity, cms to each cohort's clin table
#         results/purity_cms_scores.csv
# Time: about 10 minutes
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({
  library(tidyestimate)
  library(CMScaller)
})

run_estimate <- function(expr_symbol, is_affy) {
  df <- as.data.frame(expr_symbol)
  df <- filter_common_genes(df, id = "hgnc_symbol", tidy = FALSE, tell_missing = FALSE, find_alias = TRUE)
  out <- estimate_score(df, is_affymetrix = is_affy)
  if (!"purity" %in% colnames(out)) out$purity <- NA_real_   # RNA-seq: purity is relative only
  out
}

run_cms <- function(emat, row_names, rnaseq) {
  res <- CMScaller(emat, rowNames = row_names, RNAseq = rnaseq, doPlot = FALSE, seed = SEED)
  tibble(sample_id = rownames(res),
         cms = ifelse(is.na(res$prediction), "NC", as.character(res$prediction)))   # NC = not classified
}

all_scores <- list()

# ---- TCGA --------------------------------------------------------------------
message("\n== TCGA ==")
tcga <- readRDS(file.path(DIR_DATA, "tcga.rds"))
logtpm_sym <- collapse_by_symbol(log2(tcga$tpm + 1), tcga$genes$gene_name)
est <- run_estimate(logtpm_sym, is_affy = FALSE)
cms <- run_cms(tcga$counts, row_names = "ensg", rnaseq = TRUE)
tcga$clin <- tcga$clin %>%
  select(-any_of(c("stromal", "immune", "estimate", "purity", "cms"))) %>%
  left_join(est %>% select(sample_id = sample, stromal, immune, estimate, purity), by = "sample_id") %>%
  left_join(cms, by = "sample_id")
saveRDS(tcga, file.path(DIR_DATA, "tcga.rds"))
all_scores$TCGA <- tcga$clin
rm(tcga, logtpm_sym); gc()

# ---- GEO cohorts --------------------------------------------------------------
for (gse in GEO_COHORTS) {
  message("\n== ", gse, " ==")
  g <- readRDS(file.path(DIR_DATA, paste0(gse, ".rds")))
  est <- run_estimate(g$expr, is_affy = TRUE)
  cms <- run_cms(g$expr, row_names = "symbol", rnaseq = FALSE)
  g$clin <- g$clin %>%
    select(-any_of(c("stromal", "immune", "estimate", "purity", "cms"))) %>%
    left_join(est %>% select(sample_id = sample, stromal, immune, estimate, purity), by = "sample_id") %>%
    left_join(cms, by = "sample_id")
  saveRDS(g, file.path(DIR_DATA, paste0(gse, ".rds")))
  all_scores[[gse]] <- g$clin
}

scores <- bind_rows(lapply(all_scores, function(d)
  d %>% select(cohort, sample_id, stage4, stage_group, stromal, immune, estimate, purity, cms)))
save_table(scores, "purity_cms_scores")
print(scores %>% count(cohort, cms))
message("\n02 done.")
