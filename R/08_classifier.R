# =============================================================================
# 08_classifier.R  -  secondary: can expression separate stage I-II from III-IV? (TRIPOD)
# LASSO logistic regression on the adjusted-model genes, tuned by nested cross-validation
# in TCGA, then frozen and applied unchanged to the two GEO cohorts.
# Genes are z-scored within each cohort so the model transfers from RNA-seq to arrays.
# Output: results/classifier_performance.csv, results/classifier_genes.csv,
#         figures/figS_classifier_roc.pdf
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
})

B <- read_csv(file.path(DIR_RESULTS, "de_modelB.csv"), show_col_types = FALSE)
tc <- readRDS(file.path(DIR_DATA, "tcga_vst.rds"))
geo <- lapply(setNames(GEO_COHORTS, GEO_COHORTS), function(g) readRDS(file.path(DIR_DATA, paste0(g, ".rds"))))

# candidate features: adjusted-model genes measured on all three platforms
feat <- B %>% filter(sig, !is.na(symbol)) %>% pull(symbol) %>% unique()
feat <- Reduce(intersect, c(list(feat, rownames(tc$vst)), lapply(geo, function(x) rownames(x$expr))))
if (length(feat) < 5) {
  feat <- B %>% filter(!is.na(symbol)) %>% arrange(pvalue) %>% pull(symbol) %>% unique()
  feat <- head(Reduce(intersect, c(list(feat, rownames(tc$vst)), lapply(geo, function(x) rownames(x$expr)))), 200)
  message("Fewer than 5 significant genes shared by all platforms: using the top 200 by p-value instead.")
}
message("Features: ", length(feat))

zmat <- function(expr, ids) t(scale(t(expr[feat, ids, drop = FALSE])))   # samples x genes after transpose below
X <- t(zmat(tc$vst, tc$clin$sample_id)); X[is.na(X)] <- 0
y <- as.integer(tc$clin$stage_group == "late")

# ---- nested cross-validation in TCGA -------------------------------------------
folds <- sample(rep(seq_len(CV_FOLDS), length.out = length(y)))
cv_pred <- numeric(length(y))
for (k in seq_len(CV_FOLDS)) {
  tr <- folds != k
  m <- cv.glmnet(X[tr, ], y[tr], family = "binomial", alpha = 1, nfolds = CV_FOLDS, type.measure = "auc")
  cv_pred[!tr] <- predict(m, X[!tr, ], s = "lambda.1se", type = "response")
}
roc_cv <- roc(y, cv_pred, quiet = TRUE)

# ---- final model on all TCGA, frozen -----------------------------------------
final <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = CV_FOLDS, type.measure = "auc")
cf <- as.matrix(coef(final, s = "lambda.1se"))
genes_sel <- tibble(symbol = rownames(cf), coefficient = cf[, 1]) %>% filter(coefficient != 0, symbol != "(Intercept)") %>%
  arrange(desc(abs(coefficient)))
save_table(genes_sel, "classifier_genes")

perf <- list(tibble(cohort = "TCGA (nested CV)", n = length(y), auc = as.numeric(auc(roc_cv)),
                    ci_low = ci.auc(roc_cv)[1], ci_high = ci.auc(roc_cv)[3]))
rocs <- list(`TCGA (nested CV)` = roc_cv)
for (g in GEO_COHORTS) {
  cl <- geo[[g]]$clin
  Xg <- t(zmat(geo[[g]]$expr, cl$sample_id)); Xg[is.na(Xg)] <- 0
  yg <- as.integer(cl$stage_group == "late")
  pg <- as.numeric(predict(final, Xg, s = "lambda.1se", type = "response"))
  r <- roc(yg, pg, quiet = TRUE); rocs[[g]] <- r
  perf[[length(perf) + 1]] <- tibble(cohort = g, n = length(yg), auc = as.numeric(auc(r)),
                                     ci_low = ci.auc(r)[1], ci_high = ci.auc(r)[3])
}
perf <- bind_rows(perf) %>% mutate(genes_in_model = nrow(genes_sel))
save_table(perf, "classifier_performance")
print(perf)

p <- ggroc(rocs, legacy.axes = TRUE) + geom_abline(linetype = 2, colour = "grey60") +
  labs(colour = NULL, title = "Stage I-II vs III-IV classifier (secondary analysis)") +
  theme(legend.position = "bottom")
save_fig(p, "figS_classifier_roc", width = 6, height = 6)
message("\n08 done.")
