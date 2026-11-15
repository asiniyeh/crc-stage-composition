# =============================================================================
# 10_sensitivity.R  -  is the main finding robust? (protocol: sensitivity analyses)
# Each analysis fits an unadjusted and an adjusted limma-voom model on the SAME tumours
# and reports how many unadjusted stage genes lose significance, plus the EMT Hallmark NES.
#   S1 main specification, limma-voom instead of DESeq2 (method check)
#   S2 colon only (GSE39582 is colon only)
#   S3 published ABSOLUTE purity instead of ESTIMATE stromal/immune scores
#   S4 MSI-known tumours, MSI added to the model
#   S5 stage as an ordered trend (I-IV) instead of early/late
# Output: results/sensitivity_summary.csv
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
})

tcga <- readRDS(file.path(DIR_DATA, "tcga.rds"))
clin <- tcga$clin
if (!is.null(tcga$purity_published)) clin <- left_join(clin, tcga$purity_published, by = "sample_id")
g <- tcga$genes
keep_type <- g$gene_type %in% c("protein_coding", "lncRNA")

hall <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                 error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
pathways <- split(hall$gene_symbol, hall$gs_name)
emt <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"

fit_voom <- function(d, terms, coef) {
  cnt <- tcga$counts[keep_type, d$sample_id]
  cnt <- cnt[rowSums(cnt >= MIN_COUNT) >= ceiling(MIN_SAMPLE_FRAC * ncol(cnt)), ]
  X <- model.matrix(reformulate(terms), data = d)
  y <- calcNormFactors(DGEList(cnt))
  fit <- eBayes(lmFit(voom(y, X), X))
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  tibble(ensg = rownames(tt), symbol = g$gene_name[match(rownames(tt), g$ensg)],
         logFC = tt$logFC, t = tt$t, fdr = tt$adj.P.Val)
}
emt_nes <- function(res) {
  r <- res %>% filter(!is.na(symbol)) %>% group_by(symbol) %>% slice_max(abs(t), n = 1, with_ties = FALSE)
  fg <- fgsea(pathways[emt], setNames(r$t, r$symbol), minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE)
  c(nes = fg$NES, fdr = fg$padj)
}

run_sens <- function(label, d, termsA, termsB, coef, lfc = LFC_CUTOFF) {
  message("\n-- ", label, " (n = ", nrow(d), ")")
  A <- fit_voom(d, termsA, coef); B <- fit_voom(d, termsB, coef)
  sigA <- A$fdr < FDR_CUTOFF & abs(A$logFC) >= lfc
  sigB <- B$fdr < FDR_CUTOFF & abs(B$logFC) >= lfc
  ea <- emt_nes(A); eb <- emt_nes(B)
  tibble(analysis = label, n = nrow(d),
         n_early = if ("stage_group" %in% names(d)) sum(d$stage_group == "early") else NA,
         n_late  = if ("stage_group" %in% names(d)) sum(d$stage_group == "late") else NA,
         adjusted_for = paste(setdiff(termsB, c(termsA)), collapse = ", "),
         sig_unadjusted = sum(sigA), kept_after_adjustment = sum(sigA & sigB),
         lost_after_adjustment = sum(sigA & !sigB), pct_lost = round(100 * mean(!sigB[sigA]), 1),
         new_after_adjustment = sum(!sigA & sigB),
         EMT_NES_unadjusted = round(ea[["nes"]], 2), EMT_FDR_unadjusted = signif(ea[["fdr"]], 2),
         EMT_NES_adjusted = round(eb[["nes"]], 2), EMT_FDR_adjusted = signif(eb[["fdr"]], 2))
}

base <- clin %>%
  filter(complete.cases(stage_group, site, side, cms, stromal, immune)) %>%
  mutate(stage_group = factor(stage_group, levels = c("early", "late")),
         cms = relevel(factor(cms), ref = "CMS2"),
         stromal_z = as.numeric(scale(stromal)), immune_z = as.numeric(scale(immune)))
main_B <- c("site", "side", "cms", "stromal_z", "immune_z", "stage_group")

out <- list()
out[[1]] <- run_sens("S1 main model, limma-voom", base, c("site", "stage_group"), main_B, "stage_grouplate")

colon <- filter(base, site == "colon")
out[[2]] <- run_sens("S2 colon only", colon, "stage_group", c("side", "cms", "stromal_z", "immune_z", "stage_group"), "stage_grouplate")

if ("absolute" %in% names(base)) {
  ab <- filter(base, !is.na(absolute)) %>% mutate(absolute_z = as.numeric(scale(absolute)))
  out[[3]] <- run_sens("S3 ABSOLUTE purity instead of ESTIMATE", ab, c("site", "stage_group"),
                       c("site", "side", "cms", "absolute_z", "stage_group"), "stage_grouplate")
}

ms <- filter(base, !is.na(msi)) %>% mutate(msi = factor(msi, levels = c("MSS", "MSI-H")))
out[[4]] <- run_sens("S4 MSI-known tumours, MSI added", ms, c("site", "stage_group"),
                     c("site", "side", "msi", "cms", "stromal_z", "immune_z", "stage_group"), "stage_grouplate")

tr <- base %>% mutate(stage_num = match(stage4, c("I", "II", "III", "IV")))
# for a 1-step trend the cut-off is per stage step; 0.25 per step ~ 0.5 across early vs late
out[[5]] <- run_sens("S5 stage as ordered trend (per step)", tr, c("site", "stage_num"),
                     c("site", "side", "cms", "stromal_z", "immune_z", "stage_num"), "stage_num", lfc = LFC_CUTOFF / 2)

summ <- bind_rows(out)
save_table(summ, "sensitivity_summary")
print(as.data.frame(summ))
message("\n10 done.")
