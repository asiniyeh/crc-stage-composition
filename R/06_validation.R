# =============================================================================
# 06_validation.R  -  test the panel in GSE39582 and GSE17536 (Aim 3, Figure 5)
# Same two models as TCGA, fitted with limma on the array data:
#   A: ~ stage_group                     B: ~ covariates available + stage_group
# GSE39582 has MSI and side; GSE17536 has neither, so its model B adjusts for
# stromal, immune and CMS only (state this in the Methods).
# Output: results/validation_all_genes_<GSE>.csv, results/validation_panel.csv,
#         results/meta_analysis_panel.csv, figures/fig5_forest_panel.pdf
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({
  library(limma)
  library(metafor)
})

panel <- read_csv(file.path(DIR_RESULTS, "panel_genes.csv"), show_col_types = FALSE)
B_tcga <- read_csv(file.path(DIR_RESULTS, "de_modelB.csv"), show_col_types = FALSE)
v_tcga <- readRDS(file.path(DIR_DATA, "tcga_vst.rds"))
n_tcga <- table(v_tcga$clin$stage_group)

fit_limma <- function(expr, clin, covars) {
  vars <- c("stage_group", covars)
  cc <- complete.cases(clin[, vars])
  cl <- clin[cc, ] %>% mutate(stage_group = factor(stage_group, levels = c("early", "late")))
  if ("cms" %in% covars) cl$cms <- droplevels(relevel(factor(cl$cms), ref = "CMS2"))
  if ("stromal" %in% covars) cl$stromal <- as.numeric(scale(cl$stromal))
  if ("immune" %in% covars) cl$immune <- as.numeric(scale(cl$immune))
  f <- reformulate(c(covars, "stage_group"))
  X <- model.matrix(f, data = cl)
  fit <- eBayes(lmFit(expr[, cl$sample_id], X))
  tt <- topTable(fit, coef = "stage_grouplate", number = Inf, sort.by = "none")
  list(tab = tibble(symbol = rownames(tt), logFC = tt$logFC, t = tt$t, p = tt$P.Value, fdr = tt$adj.P.Val),
       n_early = sum(cl$stage_group == "early"), n_late = sum(cl$stage_group == "late"))
}

covars_by_cohort <- list(
  GSE39582 = c("stromal", "immune", "cms", "msi", "side"),
  GSE17536 = c("stromal", "immune", "cms")
)

val <- list(); sizes <- list(TCGA = c(early = n_tcga[["early"]], late = n_tcga[["late"]]))
for (gse in GEO_COHORTS) {
  message("\n== ", gse, " ==")
  g <- readRDS(file.path(DIR_DATA, paste0(gse, ".rds")))
  fa <- fit_limma(g$expr, g$clin, character(0))
  fb <- fit_limma(g$expr, g$clin, covars_by_cohort[[gse]])
  sizes[[gse]] <- c(early = fb$n_early, late = fb$n_late)
  all <- full_join(fa$tab %>% rename_with(~ paste0(.x, "_A"), -symbol),
                   fb$tab %>% rename_with(~ paste0(.x, "_B"), -symbol), by = "symbol")
  save_table(all, paste0("validation_all_genes_", gse))
  # does the confounding replicate? share of model-A hits lost in model B
  a_sig <- all$fdr_A < FDR_CUTOFF & abs(all$logFC_A) >= LFC_CUTOFF
  message(sprintf("  %s: %d genes significant unadjusted; %d of them stay significant after adjustment",
                  gse, sum(a_sig, na.rm = TRUE), sum(a_sig & all$fdr_B < FDR_CUTOFF, na.rm = TRUE)))
  # global agreement with TCGA model B
  m <- inner_join(B_tcga %>% filter(!is.na(symbol)) %>% distinct(symbol, .keep_all = TRUE) %>% select(symbol, statB = stat),
                  all %>% select(symbol, t_B), by = "symbol")
  message(sprintf("  Spearman correlation of adjusted statistics with TCGA: %.2f", cor(m$statB, m$t_B, method = "spearman", use = "complete")))
  val[[gse]] <- all %>% mutate(cohort = gse)
}

if (nrow(panel) == 0) {
  message("\nNo panel genes (see results/gate_result.txt): replication and meta-analysis skipped.")
} else {
# ---- panel replication -------------------------------------------------------
pv <- bind_rows(val) %>% filter(symbol %in% panel$symbol) %>%
  select(cohort, symbol, logFC_B, t_B, p_B) %>%
  left_join(panel %>% select(symbol, log2FC_tcga = log2FC), by = "symbol") %>%
  mutate(same_direction = sign(logFC_B) == sign(log2FC_tcga),
         replicated = same_direction & p_B < VALIDATION_P)
rep_tab <- pv %>% group_by(symbol) %>%
  summarise(cohorts_tested = n(), cohorts_replicated = sum(replicated, na.rm = TRUE), .groups = "drop") %>%
  mutate(validation_status = ifelse(cohorts_replicated >= 1, "replicated", "not replicated"))
save_table(pv, "validation_panel")
save_table(rep_tab, "validation_status")

# ---- random-effects meta-analysis on standardised effect sizes ---------------------
# d from a t (or Wald z) statistic: d = t * sqrt(1/n1 + 1/n2)
d_from_t <- function(t, n1, n2) t * sqrt(1 / n1 + 1 / n2)
var_d <- function(d, n1, n2) 1 / n1 + 1 / n2 + d^2 / (2 * (n1 + n2))
eff <- bind_rows(
  panel %>% transmute(symbol, cohort = "TCGA", t = stat),
  pv %>% transmute(symbol, cohort, t = t_B)
) %>% rowwise() %>%
  mutate(n1 = sizes[[cohort]][["early"]], n2 = sizes[[cohort]][["late"]],
         d = d_from_t(t, n1, n2), vi = var_d(d, n1, n2)) %>% ungroup()
meta <- eff %>% group_by(symbol) %>% filter(n() >= 2) %>%
  group_modify(~ {
    m <- rma(yi = .x$d, vi = .x$vi, method = "REML")
    tibble(pooled_d = as.numeric(m$b), ci_low = m$ci.lb, ci_high = m$ci.ub, p = m$pval, I2 = m$I2, k = m$k)
  }) %>% ungroup() %>% mutate(fdr = p.adjust(p, "BH")) %>% arrange(p)
save_table(meta, "meta_analysis_panel")

# ---- Figure 5: forest plot (top 30 genes by pooled effect) ----------------------------
top <- meta %>% slice_max(abs(pooled_d), n = 30)
fd <- bind_rows(
  eff %>% filter(symbol %in% top$symbol) %>%
    transmute(symbol, cohort, d, lo = d - 1.96 * sqrt(vi), hi = d + 1.96 * sqrt(vi)),
  top %>% transmute(symbol, cohort = "Pooled", d = pooled_d, lo = ci_low, hi = ci_high)
) %>% mutate(symbol = factor(symbol, levels = top$symbol[order(top$pooled_d)]),
             cohort = factor(cohort, levels = c("TCGA", GEO_COHORTS, "Pooled")))
p <- ggplot(fd, aes(d, symbol, colour = cohort)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.25) +
  scale_colour_manual(values = c(TCGA = "#0072B2", GSE39582 = "#E69F00", GSE17536 = "#CC79A7", Pooled = "black")) +
  labs(x = "Standardised difference, late vs early (adjusted)", y = NULL, colour = NULL,
       title = "Panel genes across three cohorts")
save_fig(p, "fig5_forest_panel", width = 7, height = max(4, 0.28 * nrow(top) + 1.5))
}
message("\n06 done.")
