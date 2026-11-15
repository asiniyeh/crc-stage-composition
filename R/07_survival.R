# =============================================================================
# 07_survival.R  -  does the panel score predict outcome beyond stage? (Aim 4, Figure 6, Table 3)
# Score (pre-specified, frozen from TCGA before validation):
#   z-score each panel gene within a cohort; score = mean z(up genes) - mean z(down genes)
# Main test: relapse-free survival in GSE39582 within stage II and within stage III.
# Output: results/table3_cox.csv, results/panel_scores.csv, figures/fig6_km_rfs.pdf
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
})

panel <- read_csv(file.path(DIR_RESULTS, "panel_genes.csv"), show_col_types = FALSE)
if (nrow(panel) == 0) stop("No panel genes (see results/gate_result.txt). Survival analysis is not run.")
up   <- panel$symbol[panel$log2FC > 0]
down <- panel$symbol[panel$log2FC < 0]

panel_score <- function(expr) {
  z <- t(scale(t(expr[intersect(c(up, down), rownames(expr)), , drop = FALSE])))
  su <- if (length(intersect(up, rownames(z))))   colMeans(z[intersect(up, rownames(z)), , drop = FALSE], na.rm = TRUE) else 0
  sd <- if (length(intersect(down, rownames(z)))) colMeans(z[intersect(down, rownames(z)), , drop = FALSE], na.rm = TRUE) else 0
  s <- su - sd
  as.numeric(scale(s))      # per 1 SD, so hazard ratios are comparable across cohorts
}

cohorts <- list(
  TCGA     = { x <- readRDS(file.path(DIR_DATA, "tcga_vst.rds")); list(expr = x$vst, clin = x$clin) },
  GSE39582 = readRDS(file.path(DIR_DATA, "GSE39582.rds")),
  GSE17536 = readRDS(file.path(DIR_DATA, "GSE17536.rds"))
)
for (nm in names(cohorts)) {
  e <- cohorts[[nm]]$expr
  cohorts[[nm]]$clin$score <- panel_score(e[, cohorts[[nm]]$clin$sample_id])
  message(sprintf("%s: %d of %d panel genes found on this platform", nm,
                  length(intersect(c(up, down), rownames(e))), length(c(up, down))))
}
save_table(bind_rows(lapply(names(cohorts), function(n)
  cohorts[[n]]$clin %>% transmute(cohort = n, sample_id, stage4, score))), "panel_scores")

fit_cox <- function(d, time, event, covars, cohort, endpoint, subset) {
  usable <- vapply(covars, function(v) v %in% names(d) && length(unique(na.omit(d[[v]]))) > 1 &&
                     mean(is.na(d[[v]])) <= 0.5, logical(1))   # same >50%-missing rule as script 04
  if (any(!usable)) message(sprintf("  %s %s %s: dropping covariates with no variation or >50%% missing: %s",
                                    cohort, endpoint, subset, paste(covars[!usable], collapse = ", ")))
  covars <- covars[usable]
  vars <- c(time, event, "score", covars)
  d <- d[complete.cases(d[, vars]), ]
  d <- d[d[[time]] > 0, ]
  if (nrow(d) < 30 || sum(d[[event]]) < 10) return(NULL)
  f  <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ", paste(c("score", covars), collapse = " + ")))
  f0 <- if (length(covars)) as.formula(paste0("Surv(", time, ", ", event, ") ~ ", paste(covars, collapse = " + "))) else NULL
  m  <- coxph(f, data = d)
  s  <- summary(m)
  c0 <- if (!is.null(f0)) summary(coxph(f0, data = d))$concordance[1] else NA
  ph <- tryCatch(cox.zph(m)$table["score", "p"], error = function(e) NA)
  tibble(cohort = cohort, endpoint = endpoint, subset = subset, adjusted_for = paste(covars, collapse = ", "),
         n = nrow(d), events = sum(d[[event]]),
         HR_per_SD = s$conf.int["score", "exp(coef)"], CI_low = s$conf.int["score", "lower .95"],
         CI_high = s$conf.int["score", "upper .95"], p = s$coefficients["score", "Pr(>|z|)"],
         c_index_without_score = c0, c_index_with_score = s$concordance[1], PH_test_p = ph)
}

res <- list()
g <- cohorts$GSE39582$clin
for (st in c("II", "III")) {
  res[[length(res) + 1]] <- fit_cox(filter(g, stage4 == st), "rfs_time", "rfs_event",
                                    c("age", "sex", "side", "msi", "chemo"), "GSE39582", "RFS", paste("stage", st))
}
res[[length(res) + 1]] <- fit_cox(filter(g, stage4 %in% c("II", "III")), "rfs_time", "rfs_event",
                                  c("stage4", "age", "sex", "side", "msi", "chemo"), "GSE39582", "RFS", "stage II-III")
res[[length(res) + 1]] <- fit_cox(g, "os_time", "os_event", c("stage4", "age", "sex", "msi"), "GSE39582", "OS", "all stages")
t1 <- cohorts$TCGA$clin
res[[length(res) + 1]] <- fit_cox(t1, "os_time", "os_event", c("stage4", "age", "sex", "msi", "site"), "TCGA", "OS", "all stages")
s1 <- cohorts$GSE17536$clin
res[[length(res) + 1]] <- fit_cox(s1, "os_time", "os_event", c("stage4", "age", "sex"), "GSE17536", "OS", "all stages")
res[[length(res) + 1]] <- fit_cox(filter(s1, stage4 %in% c("II", "III")), "rfs_time", "rfs_event",
                                  c("stage4", "age", "sex"), "GSE17536", "DFS", "stage II-III")
table3 <- bind_rows(res)
save_table(table3, "table3_cox")
print(table3)

# ---- Figure 6: Kaplan-Meier, GSE39582 RFS, median split within stage ----------
km <- list()
for (st in c("II", "III")) {
  d <- filter(g, stage4 == st, !is.na(rfs_time), !is.na(rfs_event), !is.na(score)) %>%
    mutate(group = factor(ifelse(score > median(score), "high score", "low score"), levels = c("low score", "high score")))
  fit <- survfit(Surv(rfs_time, rfs_event) ~ group, data = d)
  km[[st]] <- ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
                         palette = c("#0072B2", "#D55E00"), legend.labs = levels(d$group),
                         xlab = "Months", ylab = "Relapse-free survival",
                         title = paste0("GSE39582, stage ", st), risk.table.height = 0.25)
}
pdf(file.path(DIR_FIGS, "fig6_km_rfs.pdf"), width = 11, height = 6, onefile = FALSE)
print(arrange_ggsurvplots(km, ncol = 2, nrow = 1, print = FALSE))
dev.off()
png(file.path(DIR_FIGS, "fig6_km_rfs.png"), width = 11, height = 6, units = "in", res = 300)
print(arrange_ggsurvplots(km, ncol = 2, nrow = 1, print = FALSE))
dev.off()
message("\n07 done. The median split is for display only; Table 3 uses the continuous score.")
