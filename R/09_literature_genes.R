# =============================================================================
# 09_literature_genes.R  -  re-test the curated literature gene list (Aim 5)
# For each gene:
#   (1) tumour vs adjacent normal in TCGA (paired where possible) -> does the reported direction hold?
#   (2) late vs early stage, unadjusted (model A) and adjusted (model B)
#   (3) validation in the GEO cohorts (adjusted model)
# Status: replicated / purity-driven / not stage-associated / not measured
# Output: results/literature_genes_retest.csv
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({
  library(DESeq2)
})

lit <- read_csv(file.path(DIR_INPUT, "literature_genes.csv"), show_col_types = FALSE) %>%
  distinct(hgnc_symbol, .keep_all = TRUE)
A <- read_csv(file.path(DIR_RESULTS, "de_modelA.csv"), show_col_types = FALSE) %>% distinct(symbol, .keep_all = TRUE)
B <- read_csv(file.path(DIR_RESULTS, "de_modelB.csv"), show_col_types = FALSE) %>% distinct(symbol, .keep_all = TRUE)

# ---- (1) tumour vs normal, paired patients only ------------------------------------
tcga <- readRDS(file.path(DIR_DATA, "tcga.rds"))
nrm <- tcga$normal
pairs <- intersect(nrm$clin$patient, tcga$clin$patient)
message("Patients with tumour and adjacent normal RNA-seq: ", length(pairs))
tn <- NULL
if (length(pairs) >= 10) {
  t_ids <- tcga$clin$sample_id[match(pairs, tcga$clin$patient)]
  n_ids <- nrm$clin$sample_id[match(pairs, nrm$clin$patient)]
  cnt <- cbind(tcga$counts[, t_ids], nrm$counts[, n_ids])
  cd <- data.frame(patient = factor(c(pairs, pairs)),
                   tissue  = factor(c(rep("tumour", length(pairs)), rep("normal", length(pairs))), levels = c("normal", "tumour")))
  rownames(cd) <- colnames(cnt)
  keep <- rowSums(cnt >= MIN_COUNT) >= 3
  dds <- DESeqDataSetFromMatrix(cnt[keep, ], cd, ~ patient + tissue)
  dds <- DESeq(dds, quiet = TRUE)
  r <- results(dds, name = "tissue_tumour_vs_normal")
  tn <- tibble(ensg = rownames(r), symbol = tcga$genes$gene_name[match(rownames(r), tcga$genes$ensg)],
               tn_log2FC = r$log2FoldChange, tn_padj = r$padj) %>%
    filter(symbol %in% lit$hgnc_symbol) %>% distinct(symbol, .keep_all = TRUE)
}

# ---- (3) GEO adjusted results -------------------------------------------------------
geo <- bind_rows(lapply(GEO_COHORTS, function(g)
  read_csv(file.path(DIR_RESULTS, paste0("validation_all_genes_", g, ".csv")), show_col_types = FALSE) %>%
    transmute(symbol, cohort = g, logFC_B, p_B)))
geo_w <- geo %>% tidyr::pivot_wider(names_from = cohort, values_from = c(logFC_B, p_B))

# ---- combine ----------------------------------------------------------------
out <- lit %>%
  left_join(tn %>% select(symbol, tn_log2FC, tn_padj), by = c("hgnc_symbol" = "symbol")) %>%
  left_join(A %>% select(symbol, stageA_log2FC = log2FC, stageA_padj = padj, stageA_sig = sig), by = c("hgnc_symbol" = "symbol")) %>%
  left_join(B %>% select(symbol, stageB_log2FC = log2FC, stageB_padj = padj, stageB_sig = sig), by = c("hgnc_symbol" = "symbol")) %>%
  left_join(geo_w, by = c("hgnc_symbol" = "symbol")) %>%
  mutate(
    tumour_vs_normal = case_when(
      is.na(tn_padj) ~ "not measured",
      tn_padj >= FDR_CUTOFF ~ "no difference",
      tn_log2FC > 0 ~ "up in tumour", TRUE ~ "down in tumour"),
    reported_direction_confirmed = case_when(
      is.na(reported_direction_vs_normal) | reported_direction_vs_normal %in% c("", "NS") ~ NA,
      tumour_vs_normal == paste(tolower(reported_direction_vs_normal), "in tumour") ~ TRUE,
      TRUE ~ FALSE),
    stage_status = case_when(
      is.na(stageA_padj) & is.na(stageB_padj) ~ "not measured (filtered or absent)",
      stageB_sig %in% TRUE ~ "stage-associated after adjustment",
      stageA_sig %in% TRUE & !(stageB_sig %in% TRUE) ~ "purity-driven (lost after adjustment)",
      TRUE ~ "not stage-associated")
  ) %>% arrange(stage_status, hgnc_symbol)
save_table(out, "literature_genes_retest")
print(out %>% count(tumour_vs_normal, reported_direction_confirmed))
print(out %>% count(stage_status))
message("\n09 done.")
