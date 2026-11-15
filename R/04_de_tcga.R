# =============================================================================
# 04_de_tcga.R  -  differential expression, late vs early stage, in TCGA (Aim 1, Figure 3)
#   Model A (unadjusted): ~ site + stage_group
#   Model B (adjusted)  : ~ site + side + msi + cms + stromal + immune + stage_group
# Both models use the SAME tumours, so any difference comes from the adjustment.
# Output: results/de_modelA.csv, results/de_modelB.csv, results/de_A_vs_B_summary.csv,
#         data/tcga_vst.rds, figures/fig3_unadjusted_vs_adjusted.pdf
# Time: 20-60 minutes (DESeq2 on ~500 tumours)
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(patchwork)
  library(ggrepel)
})
workers <- max(1, parallel::detectCores() - 1)
register(if (.Platform$OS.type == "windows") SnowParam(workers) else MulticoreParam(workers))

tcga <- readRDS(file.path(DIR_DATA, "tcga.rds"))
clin <- tcga$clin

# ---- complete cases for the adjusted model ------------------------------------
model_vars <- c("stage_group", "site", "side", "msi", "cms", "stromal", "immune")
miss <- sapply(model_vars, function(v) sum(is.na(clin[[v]])))
message("Missing values per covariate:"); print(miss)
# pre-specified rule: a covariate missing in more than half of tumours is dropped, not used to drop tumours
drop_vars <- names(miss)[miss > 0.5 * nrow(clin) & names(miss) %in% c("msi", "side")]
if (length(drop_vars)) {
  warning("Dropped from model B because >50% missing: ", paste(drop_vars, collapse = ", "),
          ". Report this in the Methods.")
  model_vars <- setdiff(model_vars, drop_vars)
}
write_lines(paste("Model B covariates used:", paste(setdiff(model_vars, "stage_group"), collapse = ", ")),
            file.path(DIR_RESULTS, "de_covariates_used.txt"))
cc <- complete.cases(clin[, model_vars])
clin <- clin[cc, ] %>%
  mutate(stage_group = factor(stage_group, levels = c("early", "late")),
         site   = factor(site), side = factor(side),
         msi    = factor(msi, levels = c("MSS", "MSI-H")),
         cms    = factor(cms, levels = c("CMS2", "CMS1", "CMS3", "CMS4", "NC")),   # CMS2 = reference, the largest group
         stromal_z = as.numeric(scale(stromal)),
         immune_z  = as.numeric(scale(immune)))
clin$cms <- droplevels(clin$cms)
write_lines(sprintf("TCGA tumours in DE analysis: %d (early %d, late %d)", nrow(clin),
                    sum(clin$stage_group == "early"), sum(clin$stage_group == "late")),
            file.path(DIR_RESULTS, "de_sample_size.txt"))

# ---- genes: protein-coding + lncRNA, expressed in >= 10% of tumours -------------
cnt <- tcga$counts[, clin$sample_id]
g <- tcga$genes
keep_type <- g$gene_type %in% c("protein_coding", "lncRNA")
keep_expr <- rowSums(cnt >= MIN_COUNT) >= ceiling(MIN_SAMPLE_FRAC * ncol(cnt))
cnt <- cnt[keep_type & keep_expr, ]
g <- g[match(rownames(cnt), g$ensg), ]
message("Genes analysed: ", nrow(cnt))

run_deseq <- function(design, label) {
  message("\n-- DESeq2 ", label, ": ", deparse(design))
  cd <- as.data.frame(clin); rownames(cd) <- cd$sample_id
  stopifnot(identical(rownames(cd), colnames(cnt)))
  dds <- DESeqDataSetFromMatrix(cnt, colData = cd, design = design)
  dds <- DESeq(dds, parallel = TRUE)
  coef <- "stage_group_late_vs_early"
  res  <- results(dds, name = coef, parallel = TRUE)
  shr  <- lfcShrink(dds, coef = coef, type = "apeglm", parallel = TRUE)
  out <- tibble(
    ensg = rownames(res), symbol = g$gene_name, gene_type = g$gene_type,
    baseMean = res$baseMean, log2FC_raw = res$log2FoldChange, log2FC = shr$log2FoldChange,
    lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj
  ) %>% mutate(sig = !is.na(padj) & padj < FDR_CUTOFF & abs(log2FC) >= LFC_CUTOFF,
               direction = ifelse(sig, ifelse(log2FC > 0, "up in late", "down in late"), "ns"))
  list(dds = dds, res = out)
}

A <- run_deseq(~ site + stage_group, "model A (unadjusted)")
save_table(A$res, "de_modelA")
termsB <- c("site", intersect(c("side", "msi"), model_vars), "cms", "stromal_z", "immune_z", "stage_group")
B <- run_deseq(reformulate(termsB), "model B (adjusted)")
save_table(B$res, "de_modelB")

# variance-stabilised expression for later scripts (from the adjusted model's size factors)
vst_mat <- assay(vst(B$dds, blind = FALSE))
rownames(vst_mat) <- g$gene_name
vst_mat <- vst_mat[!duplicated(rownames(vst_mat)) & rownames(vst_mat) != "", ]
saveRDS(list(vst = vst_mat, clin = clin), file.path(DIR_DATA, "tcga_vst.rds"))

# ---- compare A and B ---------------------------------------------------------
cmp <- A$res %>% select(ensg, symbol, lfcA = log2FC, padjA = padj, sigA = sig, statA = stat) %>%
  inner_join(B$res %>% select(ensg, lfcB = log2FC, padjB = padj, sigB = sig, statB = stat), by = "ensg") %>%
  mutate(category = case_when(sigA & sigB ~ "significant in both",
                              sigA & !sigB ~ "lost after adjustment",
                              !sigA & sigB ~ "only after adjustment",
                              TRUE ~ "not significant"))
summ <- cmp %>% count(category) %>% mutate(share_of_A = ifelse(category %in% c("significant in both", "lost after adjustment"),
                                                              n / sum(cmp$sigA), NA))
save_table(summ, "de_A_vs_B_summary")
save_table(cmp, "de_A_vs_B_genes")
print(summ)

# ---- Figure 3 ----------------------------------------------------------------
volcano <- function(res, title) {
  top <- res %>% filter(sig) %>% arrange(padj) %>% head(12)
  ggplot(res, aes(log2FC, -log10(pvalue), colour = direction)) +
    geom_point(size = 0.5, alpha = 0.6) +
    geom_text_repel(data = top, aes(label = symbol), size = 2.6, colour = "black", max.overlaps = 20) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = 2, colour = "grey50") +
    scale_colour_manual(values = c("up in late" = "#D55E00", "down in late" = "#0072B2", ns = "grey75")) +
    labs(x = "log2 fold change, late vs early (shrunk)", y = "-log10 p", colour = NULL, title = title) +
    theme(legend.position = "bottom")
}
pc <- ggplot(cmp, aes(lfcA, lfcB, colour = category)) +
  geom_abline(linetype = 2, colour = "grey50") + geom_hline(yintercept = 0, colour = "grey80") +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_colour_manual(values = c("significant in both" = "#009E73", "lost after adjustment" = "#D55E00",
                                 "only after adjustment" = "#0072B2", "not significant" = "grey80")) +
  labs(x = "log2 FC, model A (unadjusted)", y = "log2 FC, model B (adjusted)", colour = NULL,
       title = paste0("C  Effect of adjusting for ", paste(setdiff(model_vars, "stage_group"), collapse = ", "))) +
  theme(legend.position = "bottom")
save_fig((volcano(A$res, "A  Unadjusted") | volcano(B$res, "B  Adjusted")) / pc,
         "fig3_unadjusted_vs_adjusted", width = 11, height = 10)
message("\n04 done.")
