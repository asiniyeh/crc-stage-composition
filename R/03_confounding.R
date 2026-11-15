# =============================================================================
# 03_confounding.R  -  does tumour composition change with stage?  (Aim 1, Figure 2)
# Output: results/confounding_tests.csv, figures/fig2_composition_by_stage.pdf
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages(library(patchwork))

scores <- read_csv(file.path(DIR_RESULTS, "purity_cms_scores.csv"), show_col_types = FALSE) %>%
  filter(!is.na(stage4)) %>%
  mutate(stage4 = factor(stage4, levels = c("I", "II", "III", "IV")),
         stage_num = as.integer(stage4),
         cohort = factor(cohort, levels = c("TCGA", GEO_COHORTS)))

tests <- list()
for (co in levels(scores$cohort)) {
  d <- filter(scores, cohort == co)
  for (v in c("stromal", "immune", "estimate")) {
    kw <- kruskal.test(d[[v]] ~ d$stage4)
    sp <- suppressWarnings(cor.test(d[[v]], d$stage_num, method = "spearman"))
    tests[[length(tests) + 1]] <- tibble(
      cohort = co, variable = v, test = "Kruskal-Wallis across stage I-IV",
      statistic = unname(kw$statistic), p = kw$p.value)
    tests[[length(tests) + 1]] <- tibble(
      cohort = co, variable = v, test = "Spearman trend with stage",
      statistic = unname(sp$estimate), p = sp$p.value)
  }
  tab <- table(d$cms, d$stage_group)
  ch <- suppressWarnings(chisq.test(tab))
  tests[[length(tests) + 1]] <- tibble(
    cohort = co, variable = "CMS", test = "Chi-square CMS x early/late",
    statistic = unname(ch$statistic), p = ch$p.value)
  # proportion CMS4 early vs late, the key number for the Results text
  p4 <- d %>% group_by(stage_group) %>% summarise(cms4 = mean(cms == "CMS4"), .groups = "drop")
  tests[[length(tests) + 1]] <- tibble(
    cohort = co, variable = "CMS4 share",
    test = paste0("early ", round(100 * p4$cms4[p4$stage_group == "early"]), "% vs late ",
                  round(100 * p4$cms4[p4$stage_group == "late"]), "%"),
    statistic = NA_real_, p = NA_real_)
}
tests <- bind_rows(tests) %>% mutate(p_fdr = p.adjust(p, "BH"))
save_table(tests, "confounding_tests")
print(tests, n = Inf)

# ---- Figure 2 ----------------------------------------------------------------
pa <- ggplot(scores, aes(stage4, stromal)) +
  geom_boxplot(outlier.size = 0.6, fill = "grey90") +
  facet_wrap(~ cohort, scales = "free_y") +
  labs(x = "AJCC stage", y = "ESTIMATE stromal score", title = "A  Stromal content by stage")
pb <- ggplot(scores, aes(stage4, immune)) +
  geom_boxplot(outlier.size = 0.6, fill = "grey90") +
  facet_wrap(~ cohort, scales = "free_y") +
  labs(x = "AJCC stage", y = "ESTIMATE immune score", title = "B  Immune content by stage")
pc <- scores %>% count(cohort, stage4, cms) %>% group_by(cohort, stage4) %>% mutate(share = n / sum(n)) %>%
  ggplot(aes(stage4, share, fill = cms)) +
  geom_col() + facet_wrap(~ cohort) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(CMS1 = "#E69F00", CMS2 = "#0072B2", CMS3 = "#CC79A7", CMS4 = "#009E73", NC = "grey75")) +
  labs(x = "AJCC stage", y = "Share of tumours", fill = "CMS", title = "C  Consensus molecular subtype by stage")
save_fig(pa / pb / pc, "fig2_composition_by_stage", width = 9, height = 11)
message("\n03 done.")
