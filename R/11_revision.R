# =============================================================================
# 11_revision.R  -  analyses added after internal pre-submission review
#   R1 attenuation of effect sizes (DESeq2 models A vs B), not just loss of significance
#   R2 collinearity: how much adjustment inflates the SE of the stage coefficient
#   R3 decomposition of the adjustment (side / scores only / CMS only / both), limma-voom
#   R4 null benchmark: adjust for random covariates of the same dimension
#   R5 classifier with gene selection inside the cross-validation loop
#   R6 GSE39582 with harmonised covariates (no MSI, no side)
#   R7 permutation null for the panel replication count; Hartung-Knapp meta-analysis
#   R8 tumours excluded by complete-case analysis
# Output: results/rev_*.csv
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({
  library(limma); library(edgeR); library(fgsea); library(msigdbr)
  library(glmnet); library(pROC); library(metafor)
})
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)

# ---------------------------------------------------------------- R1 attenuation
A <- read_csv(file.path(DIR_RESULTS, "de_modelA.csv"), show_col_types = FALSE)
B <- read_csv(file.path(DIR_RESULTS, "de_modelB.csv"), show_col_types = FALSE)
m <- inner_join(A, B, by = "ensg", suffix = c("A", "B")) %>% filter(sigA)
cat_lab <- ifelse(m$sigB, "kept", "lost")
att <- function(d, lab) {
  ratio <- d$log2FCB / d$log2FCA
  slope <- coef(lm(log2FCB ~ 0 + log2FCA, data = d))[[1]]
  excl  <- abs(d$log2FC_rawA - d$log2FC_rawB) > 1.96 * d$lfcSEB
  tibble(set = lab, n = nrow(d), median_ratio_B_over_A = median(ratio), slope_B_on_A = slope,
         pct_shrunk_over_half = 100 * mean(ratio < 0.5), pct_sign_flip = 100 * mean(ratio < 0),
         pct_A_estimate_outside_B_95CI = 100 * mean(excl),
         median_SE_ratio_B_over_A = median(d$lfcSEB / d$lfcSEA))
}
r1 <- bind_rows(att(m, "all 760 unadjusted stage genes"),
                att(m[cat_lab == "lost", ], "lost after adjustment"),
                att(m[cat_lab == "kept", ], "kept after adjustment"))
save_table(r1, "rev_R1_attenuation"); print(as.data.frame(r1))

# ---------------------------------------------------------------- data for limma analyses
tcga <- readRDS(file.path(DIR_DATA, "tcga.rds"))
g <- tcga$genes
clin_all <- tcga$clin
base <- clin_all %>%
  filter(complete.cases(stage_group, site, side, cms, stromal, immune)) %>%
  mutate(stage_group = factor(stage_group, levels = c("early", "late")),
         cms = relevel(factor(cms), ref = "CMS2"),
         stromal_z = as.numeric(scale(stromal)), immune_z = as.numeric(scale(immune)))
cnt <- tcga$counts[g$gene_type %in% c("protein_coding", "lncRNA"), base$sample_id]
cnt <- cnt[rowSums(cnt >= MIN_COUNT) >= ceiling(MIN_SAMPLE_FRAC * ncol(cnt)), ]
dge <- calcNormFactors(DGEList(cnt))
sym <- g$gene_name[match(rownames(cnt), g$ensg)]
full_terms <- c("site", "side", "cms", "stromal_z", "immune_z", "stage_group")
Xfull <- model.matrix(reformulate(full_terms), base)
v <- voom(dge, Xfull)                      # weights from the full model, reused below (approximation for R4)

hall <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                 error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
emt <- list(EMT = hall$gene_symbol[hall$gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"])

fit_terms <- function(terms, d = base, vv = NULL) {
  X <- model.matrix(reformulate(terms), d)
  if (is.null(vv)) vv <- voom(dge, X)
  fit <- eBayes(lmFit(vv, X))
  tt <- topTable(fit, coef = "stage_grouplate", number = Inf, sort.by = "none")
  tibble(ensg = rownames(tt), symbol = sym, logFC = tt$logFC, t = tt$t, fdr = tt$adj.P.Val,
         sig = tt$adj.P.Val < FDR_CUTOFF & abs(tt$logFC) >= LFC_CUTOFF)
}
emt_nes <- function(r) {
  r <- r %>% filter(!is.na(symbol)) %>% group_by(symbol) %>% slice_max(abs(t), n = 1, with_ties = FALSE)
  f <- fgsea(emt, setNames(r$t, r$symbol), minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE)
  c(nes = f$NES, fdr = f$padj)
}
summ_vs <- function(ref, r, lab) {
  e <- emt_nes(r); s <- ref$sig
  tibble(model = lab, sig_genes = sum(r$sig), ref_genes_lost = sum(s & !r$sig),
         pct_ref_lost = round(100 * mean(!r$sig[s]), 1),
         slope_vs_ref = round(coef(lm(r$logFC[s] ~ 0 + ref$logFC[s]))[[1]], 2),
         EMT_NES = round(e[["nes"]], 2), EMT_FDR = signif(e[["fdr"]], 2))
}

# ---------------------------------------------------------------- R2 collinearity
cov_only <- model.matrix(reformulate(setdiff(full_terms, "stage_group")), base)[, -1]
r2_stage <- summary(lm(as.numeric(base$stage_group == "late") ~ cov_only))$r.squared
r2_siteonly <- summary(lm(as.numeric(base$stage_group == "late") ~ base$site))$r.squared
cms4 <- as.numeric(base$cms == "CMS4")
r2 <- tibble(
  quantity = c("R2 of stage on model-B covariates", "VIF of stage in model B",
               "Expected SE inflation for stage, model B vs model A",
               "Point-biserial r, stromal score vs CMS4", "R2 of stromal score on CMS"),
  value = c(r2_stage, 1 / (1 - r2_stage), sqrt((1 - r2_siteonly) / (1 - r2_stage)),
            cor(base$stromal_z, cms4), summary(lm(stromal_z ~ cms, base))$r.squared))
save_table(r2, "rev_R2_collinearity"); print(as.data.frame(r2))

# ---------------------------------------------------------------- R3 decomposition
message("R3 decomposition")
mods <- list(
  "A  site (unadjusted)"            = c("site", "stage_group"),
  "A1 + side"                       = c("site", "side", "stage_group"),
  "C  + side + stromal + immune"    = c("site", "side", "stromal_z", "immune_z", "stage_group"),
  "D  + side + CMS"                 = c("site", "side", "cms", "stage_group"),
  "B  + side + CMS + stromal + immune" = full_terms)
fits <- lapply(mods, fit_terms)
r3 <- bind_rows(lapply(names(fits), function(n) summ_vs(fits[[1]], fits[[n]], n)))
save_table(r3, "rev_R3_decomposition"); print(as.data.frame(r3))

# ---------------------------------------------------------------- R4 null benchmark
message("R4 null benchmark")
ref <- fit_terms(mods[[2]], vv = v)        # site + side, same voom weights as the null fits
obsB <- fit_terms(full_terms, vv = v)
k <- ncol(Xfull) - ncol(model.matrix(reformulate(mods[[2]]), base))   # df added by CMS + scores
vst_like <- v$E
nullrun <- function(type) {
  d <- base
  Z <- if (type == "gaussian") matrix(rnorm(nrow(d) * k), ncol = k) else
    t(scale(t(vst_like[sample(nrow(vst_like), k), , drop = FALSE])))  |> t()
  colnames(Z) <- paste0("z", seq_len(k)); d <- bind_cols(d, as_tibble(Z))
  r <- fit_terms(c("site", "side", colnames(Z), "stage_group"), d = d, vv = v)
  c(pct_lost = 100 * mean(!r$sig[ref$sig]), slope = coef(lm(r$logFC[ref$sig] ~ 0 + ref$logFC[ref$sig]))[[1]])
}
NPERM <- 50   # weighted lmFit is slow on 497 tumours; 50 per null type
nulls <- bind_rows(
  bind_rows(parallel::mclapply(seq_len(NPERM), function(i) as_tibble_row(nullrun("gaussian")), mc.cores = 2, mc.set.seed = TRUE)) %>% mutate(type = "random Gaussian covariates"),
  bind_rows(parallel::mclapply(seq_len(NPERM), function(i) as_tibble_row(nullrun("genes")), mc.cores = 2, mc.set.seed = TRUE)) %>% mutate(type = "expression of random genes"))
obs <- c(pct_lost = 100 * mean(!obsB$sig[ref$sig]), slope = coef(lm(obsB$logFC[ref$sig] ~ 0 + ref$logFC[ref$sig]))[[1]])
r4 <- nulls %>% group_by(type) %>%
  summarise(k_covariates = k, null_pct_lost_median = median(pct_lost),
            null_pct_lost_95th = quantile(pct_lost, 0.95), null_pct_lost_max = max(pct_lost),
            null_slope_median = median(slope), null_slope_5th = quantile(slope, 0.05),
            observed_pct_lost = obs[["pct_lost"]], observed_slope = obs[["slope"]],
            perm_p_pct_lost = (sum(pct_lost >= obs[["pct_lost"]]) + 1) / (n() + 1), .groups = "drop")
save_table(r4, "rev_R4_null_benchmark"); print(as.data.frame(r4))

# ---------------------------------------------------------------- R5 classifier, selection inside CV
message("R5 classifier")
tc <- readRDS(file.path(DIR_DATA, "tcga_vst.rds"))
geo <- lapply(setNames(GEO_COHORTS, GEO_COHORTS), function(x) readRDS(file.path(DIR_DATA, paste0(x, ".rds"))))
shared <- Reduce(intersect, c(list(rownames(tc$vst)), lapply(geo, function(x) rownames(x$expr))))
cl <- tc$clin; stopifnot(all(cl$sample_id %in% base$sample_id))
idx <- match(cl$sample_id, base$sample_id)
y <- as.integer(cl$stage_group == "late")
zs <- function(M) { Z <- t(scale(t(M))); Z[is.na(Z)] <- 0; t(Z) }
select_genes <- function(rows) {
  d <- base[idx[rows], ]
  X <- model.matrix(reformulate(full_terms), d); X <- X[, colSums(abs(X)) > 0, drop = FALSE]
  dd <- calcNormFactors(DGEList(cnt[, d$sample_id]))
  fit <- eBayes(lmFit(voom(dd, X), X))
  tt <- topTable(fit, coef = "stage_grouplate", number = Inf, sort.by = "none")
  s <- unique(sym[tt$adj.P.Val < FDR_CUTOFF & abs(tt$logFC) >= LFC_CUTOFF])
  intersect(s, shared)
}
folds <- sample(rep(seq_len(CV_FOLDS), length.out = length(y)))
pred <- numeric(length(y)); nfeat <- integer(CV_FOLDS)
for (kk in seq_len(CV_FOLDS)) {
  tr <- folds != kk
  f <- select_genes(which(tr)); nfeat[kk] <- length(f)
  Xtr <- zs(tc$vst[f, cl$sample_id[tr]]); Xte <- zs(tc$vst[f, cl$sample_id[!tr]])
  mm <- cv.glmnet(Xtr, y[tr], family = "binomial", alpha = 1, nfolds = CV_FOLDS, type.measure = "auc")
  pred[!tr] <- predict(mm, Xte, s = "lambda.1se", type = "response")
}
rc <- roc(y, pred, quiet = TRUE)
calib <- function(yy, p) { lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  c(slope = coef(glm(yy ~ lp, family = binomial))[[2]],
    citl = coef(glm(yy ~ offset(lp), family = binomial))[[1]]) }
# frozen model from 08 (features selected on all TCGA) for external calibration
feat <- B %>% filter(sig, !is.na(symbol)) %>% pull(symbol) %>% unique() %>% intersect(shared)
final <- cv.glmnet(zs(tc$vst[feat, cl$sample_id]), y, family = "binomial", alpha = 1, nfolds = CV_FOLDS, type.measure = "auc")
r5 <- list(tibble(cohort = "TCGA, selection inside CV", n = length(y), auc = as.numeric(auc(rc)),
                  ci_low = ci.auc(rc)[1], ci_high = ci.auc(rc)[3],
                  calib_slope = calib(y, pred)[["slope"]], calib_in_large = calib(y, pred)[["citl"]],
                  features_per_fold = paste(nfeat, collapse = "/")))
for (gg in GEO_COHORTS) {
  c2 <- geo[[gg]]$clin; yg <- as.integer(c2$stage_group == "late")
  pg <- as.numeric(predict(final, zs(geo[[gg]]$expr[feat, c2$sample_id]), s = "lambda.1se", type = "response"))
  rr <- roc(yg, pg, quiet = TRUE); cb <- calib(yg, pg)
  r5[[length(r5) + 1]] <- tibble(cohort = paste(gg, "(frozen model)"), n = length(yg), auc = as.numeric(auc(rr)),
                                 ci_low = ci.auc(rr)[1], ci_high = ci.auc(rr)[3],
                                 calib_slope = cb[["slope"]], calib_in_large = cb[["citl"]], features_per_fold = NA)
}
r5 <- bind_rows(r5); save_table(r5, "rev_R5_classifier"); print(as.data.frame(r5))

# ---------------------------------------------------------------- R6 + R7 validation
message("R6/R7 validation")
panel <- read_csv(file.path(DIR_RESULTS, "panel_genes.csv"), show_col_types = FALSE)
fit_limma <- function(expr, cl, covars, perm = FALSE) {
  cc <- complete.cases(cl[, c("stage_group", covars)]); cl <- cl[cc, ]
  cl$stage_group <- factor(cl$stage_group, levels = c("early", "late"))
  if (perm) cl$stage_group <- sample(cl$stage_group)
  if ("cms" %in% covars) cl$cms <- droplevels(relevel(factor(cl$cms), ref = "CMS2"))
  for (s in intersect(c("stromal", "immune"), covars)) cl[[s]] <- as.numeric(scale(cl[[s]]))
  X <- model.matrix(reformulate(c(covars, "stage_group")), cl)
  tt <- topTable(eBayes(lmFit(expr[, cl$sample_id], X)), coef = "stage_grouplate", number = Inf, sort.by = "none")
  tibble(symbol = rownames(tt), logFC = tt$logFC, t = tt$t, p = tt$P.Value)
}
cov_main <- list(GSE39582 = c("stromal", "immune", "cms", "msi", "side"), GSE17536 = c("stromal", "immune", "cms"))
cov_harm <- list(GSE39582 = c("stromal", "immune", "cms"), GSE17536 = c("stromal", "immune", "cms"))
pdir <- setNames(sign(panel$log2FC), panel$symbol)
count_rep <- function(res) {
  rep <- lapply(res, function(r) { r <- r[r$symbol %in% names(pdir), ]
    setNames(sign(r$logFC) == pdir[r$symbol] & r$p < VALIDATION_P, r$symbol) })
  gs <- intersect(names(rep[[1]]), names(rep[[2]]))
  c(tested = length(gs), any = sum(rep[[1]][gs] | rep[[2]][gs]), both = sum(rep[[1]][gs] & rep[[2]][gs]))
}
obs_main <- count_rep(lapply(GEO_COHORTS, function(x) fit_limma(geo[[x]]$expr, geo[[x]]$clin, cov_main[[x]])))
obs_harm <- count_rep(lapply(GEO_COHORTS, function(x) fit_limma(geo[[x]]$expr, geo[[x]]$clin, cov_harm[[x]])))
NP <- 500
perm_counts <- t(replicate(NP, count_rep(lapply(GEO_COHORTS, function(x)
  fit_limma(geo[[x]]$expr, geo[[x]]$clin, cov_main[[x]], perm = TRUE)))))
analytic <- obs_main[["tested"]] * (1 - (1 - VALIDATION_P / 2)^2)
r67 <- tibble(
  analysis = c("Main covariates (GSE39582 incl. MSI, side)", "Harmonised covariates (stromal, immune, CMS)",
               "Permutation null, main covariates (median [95th pct])", "Analytic expectation (independence)"),
  genes_tested = c(obs_main[["tested"]], obs_harm[["tested"]], obs_main[["tested"]], obs_main[["tested"]]),
  replicated_in_any = c(obs_main[["any"]], obs_harm[["any"]],
                        sprintf("%g [%g]", median(perm_counts[, "any"]), quantile(perm_counts[, "any"], 0.95)),
                        sprintf("%.1f", analytic)),
  replicated_in_both = c(obs_main[["both"]], obs_harm[["both"]],
                         sprintf("%g [%g]", median(perm_counts[, "both"]), quantile(perm_counts[, "both"], 0.95)), "0.1"),
  perm_p_any = c((sum(perm_counts[, "any"] >= obs_main[["any"]]) + 1) / (NP + 1), NA, NA, NA))
save_table(r67, "rev_R6_R7_replication"); print(as.data.frame(r67))

meta <- read_csv(file.path(DIR_RESULTS, "meta_analysis_panel.csv"), show_col_types = FALSE)
eff <- read_csv(file.path(DIR_RESULTS, "validation_panel.csv"), show_col_types = FALSE)
# rebuild per-cohort d as in 06 and refit with Hartung-Knapp
v_tcga <- readRDS(file.path(DIR_DATA, "tcga_vst.rds")); nt <- table(v_tcga$clin$stage_group)
sizes <- list(TCGA = c(nt[["early"]], nt[["late"]]))
for (x in GEO_COHORTS) { c2 <- geo[[x]]$clin[complete.cases(geo[[x]]$clin[, c("stage_group", cov_main[[x]])]), ]
  sizes[[x]] <- c(sum(c2$stage_group == "early"), sum(c2$stage_group == "late")) }
E <- bind_rows(panel %>% transmute(symbol, cohort = "TCGA", t = stat), eff %>% transmute(symbol, cohort, t = t_B)) %>%
  rowwise() %>% mutate(n1 = sizes[[cohort]][1], n2 = sizes[[cohort]][2], d = t * sqrt(1 / n1 + 1 / n2),
                       vi = 1 / n1 + 1 / n2 + d^2 / (2 * (n1 + n2))) %>% ungroup()
hk <- E %>% group_by(symbol) %>% filter(n() >= 3) %>% group_modify(~ {
  mm <- rma(yi = .x$d, vi = .x$vi, method = "REML", test = "knha"); tibble(d = as.numeric(mm$b), p_hk = mm$pval) }) %>%
  ungroup() %>% mutate(fdr_hk = p.adjust(p_hk, "BH"))
r7b <- tibble(method = c("REML, Wald (as in manuscript)", "REML, Hartung-Knapp"),
              genes_meta = c(sum(meta$k >= 3), nrow(hk)),
              fdr_below_0.05 = c(sum(meta$fdr < 0.05 & meta$k >= 3), sum(hk$fdr_hk < 0.05)),
              genes = c(paste(meta$symbol[meta$fdr < 0.05 & meta$k >= 3], collapse = ", "),
                        paste(hk$symbol[hk$fdr_hk < 0.05], collapse = ", ")))
save_table(r7b, "rev_R7_meta_hk"); print(as.data.frame(r7b))

# ---------------------------------------------------------------- R8 excluded tumours
ex <- clin_all %>% mutate(included = sample_id %in% base$sample_id)
r8 <- bind_rows(
  ex %>% count(included, stage4) %>% mutate(variable = "stage", level = stage4) %>% select(-stage4),
  ex %>% count(included, cms) %>% mutate(variable = "CMS", level = cms) %>% select(-cms),
  ex %>% count(included, site) %>% mutate(variable = "site", level = site) %>% select(-site),
  ex %>% filter(!included) %>% summarise(n = sum(is.na(side)), included = FALSE) %>% mutate(variable = "missing", level = "side"),
  ex %>% filter(!included) %>% summarise(n = sum(is.na(cms)), included = FALSE) %>% mutate(variable = "missing", level = "cms"),
  ex %>% filter(!included) %>% summarise(n = sum(is.na(stromal)), included = FALSE) %>% mutate(variable = "missing", level = "stromal")) %>%
  group_by(variable, included) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
r8_test <- tibble(test = c("stage, chi-sq", "CMS, chi-sq"),
  p = c(chisq.test(table(ex$included, ex$stage4))$p.value, chisq.test(table(ex$included, ex$cms))$p.value))
save_table(r8, "rev_R8_excluded"); save_table(r8_test, "rev_R8_excluded_tests")
print(as.data.frame(r8)); print(r8_test)
message("\n11 done.")
