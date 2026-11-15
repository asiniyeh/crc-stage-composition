# =============================================================================
# 05_panel_gsea.R  -  epithelium-intrinsic panel (Aim 2) and Hallmark GSEA (Figure 4)
# Output: results/panel_genes.csv, results/gsea_hallmark.csv,
#         figures/fig4_hallmark_gsea.pdf
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(tidyestimate)
})

A <- read_csv(file.path(DIR_RESULTS, "de_modelA.csv"), show_col_types = FALSE)
B <- read_csv(file.path(DIR_RESULTS, "de_modelB.csv"), show_col_types = FALSE)
v <- readRDS(file.path(DIR_DATA, "tcga_vst.rds"))

# ---- panel -------------------------------------------------------------------
estimate_genes <- unique(c(tidyestimate::gene_sets$stromal_signature, tidyestimate::gene_sets$immune_signature))
cand <- B %>% filter(sig, !is.na(symbol), symbol %in% rownames(v$vst))
rho <- vapply(cand$symbol, function(s)
  suppressWarnings(cor(v$vst[s, ], v$clin$stromal, method = "spearman")), numeric(1))
cand <- cand %>%
  mutate(rho_stromal = rho,
         in_estimate_sets = symbol %in% estimate_genes,
         panel = !in_estimate_sets & abs(rho_stromal) < STROMA_RHO_MAX) %>%
  left_join(A %>% select(ensg, log2FC_A = log2FC, padj_A = padj), by = "ensg") %>%
  arrange(desc(panel), padj)
save_table(cand, "panel_candidates")
panel <- filter(cand, panel)
save_table(panel, "panel_genes")

gate <- sprintf("GATE: %d genes significant in model B; %d pass the stroma filter and form the panel.",
                nrow(cand), nrow(panel))
if (nrow(panel) == 0) gate <- paste(gate, "No panel: frame the paper around stromal confounding (see protocol).")
write_lines(gate, file.path(DIR_RESULTS, "gate_result.txt"))
message("\n", gate)

# ---- Hallmark GSEA, both models ---------------------------------------------------
hall <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                 error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
pathways <- split(hall$gene_symbol, hall$gs_name)

run_gsea <- function(res, label) {
  r <- res %>% filter(!is.na(stat), !is.na(symbol), symbol != "") %>%
    group_by(symbol) %>% slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
  stats <- setNames(r$stat, r$symbol)
  fg <- fgsea(pathways, stats, minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE, eps = 0)
  as_tibble(fg) %>% mutate(model = label, leadingEdge = vapply(leadingEdge, paste, "", collapse = ";"))
}
gs <- bind_rows(run_gsea(A, "A unadjusted"), run_gsea(B, "B adjusted"))
save_table(gs %>% select(model, pathway, NES, pval, padj, size, leadingEdge), "gsea_hallmark")

keep_paths <- gs %>% group_by(pathway) %>% filter(any(padj < FDR_CUTOFF)) %>% pull(pathway) %>% unique()
pd <- gs %>% filter(pathway %in% keep_paths) %>%
  mutate(pathway = gsub("_", " ", sub("^HALLMARK_", "", pathway)),
         pathway = reorder(pathway, NES * (model == "A unadjusted")),
         sig = padj < FDR_CUTOFF)
p <- ggplot(pd, aes(model, pathway)) +
  geom_point(aes(size = -log10(padj), colour = NES, shape = sig)) +
  scale_colour_gradient2(low = "#0072B2", mid = "grey90", high = "#D55E00") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "not significant")) +
  labs(x = NULL, y = NULL, shape = NULL, size = "-log10 FDR",
       title = "Hallmark pathways, late vs early stage")
save_fig(p, "fig4_hallmark_gsea", width = 7.5, height = max(4, 0.25 * length(keep_paths) + 2))
message("\n05 done.")
