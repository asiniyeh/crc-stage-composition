# =============================================================================
# run_all.R  -  run the whole pipeline in order
# Open R or RStudio IN THIS FOLDER, then:   source("run_all.R")
# If a step fails, fix it and re-run from that step: each script saves its own outputs.
# =============================================================================
scripts <- c(
  "R/01_download_clean.R",
  "R/02_purity_cms.R",
  "R/03_confounding.R",
  "R/04_de_tcga.R",
  "R/05_panel_gsea.R",
  "R/06_validation.R",
  "R/07_survival.R",
  "R/08_classifier.R",
  "R/09_literature_genes.R"
)
start_at <- 1          # change to e.g. 4 to resume from script 04

dir.create("logs", showWarnings = FALSE)
for (s in scripts[start_at:length(scripts)]) {
  message("\n==================== ", s, " ====================")
  log <- file.path("logs", sub("\\.R$", ".log", basename(s)))
  t0 <- Sys.time()
  ok <- tryCatch({
    withCallingHandlers(
      source(s, echo = FALSE, local = new.env()),
      message = function(m) cat(conditionMessage(m), file = log, append = TRUE),
      warning = function(w) cat("WARNING: ", conditionMessage(w), "\n", file = log, append = TRUE)
    )
    TRUE
  }, error = function(e) {
    cat("ERROR: ", conditionMessage(e), "\n", file = log, append = TRUE)
    message("\nSTOPPED in ", s, ":\n  ", conditionMessage(e),
            "\nSend logs/", basename(log), " to Claude.")
    FALSE
  })
  message(sprintf("  %s finished in %.1f min", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (!ok) break
}
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
message("\nDone. Zip the 'results' and 'figures' folders (and 'logs') and send them back.")
