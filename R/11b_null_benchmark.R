# =============================================================================
# 11_revision.R  -  analyses added after internal pre-submission review
# 11b: R4 null benchmark with 200 draws per type (re-run of R4 in 11_revision.R)
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

mods <- list(c("site","stage_group"), c("site","side","stage_group"))

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
NPERM <- 100
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
save_table(r4, "rev_R4_null_benchmark"); save_table(nulls, "rev_R4_null_draws"); print(as.data.frame(r4))


message("11b done.")
