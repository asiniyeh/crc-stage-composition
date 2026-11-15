# Manuscript 1 – Stage-associated gene expression in colorectal cancer, adjusted for tumour purity

Analysis pipeline for the study protocol "Manuscript 1 – Stage-associated transcriptomic changes in colorectal cancer".
Everything runs on public data: TCGA-COAD/READ, GSE39582 and GSE17536.

## What you need

- **R 4.3 or newer** (https://cran.r-project.org) and **RStudio** (optional but easier)
- **16 GB RAM** recommended (8 GB may work but DESeq2 will be slow)
- **About 6 GB free disk space** (TCGA downloads are 2–3 GB)
- **Internet access** for the first run (TCGA, GEO, CRAN, Bioconductor, GitHub)

## How to run

1. Unzip this folder somewhere simple, e.g. `Documents/m1_pipeline`.
2. Open RStudio → *File → Open Project in New Session*… or simply set the working directory to this folder:
   ```r
   setwd("~/Documents/m1_pipeline")
   ```
3. Install packages once (20–40 min):
   ```r
   source("R/00_install.R")
   ```
4. Run the whole analysis (2–4 hours, mostly downloading and DESeq2):
   ```r
   source("run_all.R")
   ```
   If it stops, the console says which script failed and the log file to send. After fixing, set
   `start_at` in `run_all.R` to that script's number and run it again; earlier steps are not repeated.
5. Zip the folders **results**, **figures** and **logs** and send them back. Do not send `data_raw/` (too large).

## The scripts

| Script | What it does | Main outputs |
| --- | --- | --- |
| 01_download_clean | Downloads TCGA + GEO, harmonises clinical data | Table 1, flow counts |
| 02_purity_cms | ESTIMATE stromal/immune scores, CMS calls | purity_cms_scores.csv |
| 03_confounding | Composition vs stage | Figure 2, confounding_tests.csv |
| 04_de_tcga | DESeq2, unadjusted vs adjusted | Figure 3, de_modelA/B.csv |
| 05_panel_gsea | Epithelium-intrinsic panel, Hallmark GSEA, **gate result** | Figure 4, panel_genes.csv, gate_result.txt |
| 06_validation | limma in the GEO cohorts, meta-analysis | Figure 5, validation_*.csv |
| 07_survival | Panel score vs RFS/OS | Figure 6, Table 3 |
| 08_classifier | LASSO stage classifier (secondary) | ROC figure, classifier_performance.csv |
| 09_literature_genes | Re-tests the curated literature gene list | literature_genes_retest.csv |
| 10_sensitivity | Five limma-voom sensitivity analyses (Table 2) | sensitivity_summary.csv |
| 11_revision | Effect-size attenuation, collinearity, decomposition of the adjustment, classifier with selection inside CV, replication permutation null, Hartung-Knapp meta-analysis, excluded tumours (Supplementary Table S6) | rev_R*.csv |
| 11b_null_benchmark | Null benchmark with random covariates, 100 draws per type (values reported) | rev_R4_null_benchmark.csv, rev_R4_null_draws.csv |

All thresholds are fixed in `R/00_config.R` and match the protocol.

## Notes and known issues

- **GEO column names.** Script 01 finds clinical fields by pattern (e.g. `tnm.stage`, `rfs.event`). If GEO has renamed one, the script stops and prints every available column; change the pattern on the line it names and re-run from script 01.
- **Event coding.** At the end of script 01, check the printed event counts (0 = censored, 1 = event) look sensible before trusting survival results.
- **GSE17536** has no MSI or tumour-side data, so its adjusted model uses stromal, immune and CMS only.
- **Literature gene list** (`input/literature_genes.csv`): old names were updated to current HGNC symbols (IL8→CXCL8, L-FABP→FABP1, CXCR7→ACKR3, KIAA1199→CEMIP, and others). Two Nosho et al. mappings are probable rather than certain (Rho 8→RND3, P120→CTNND1). Entries that could not be mapped were left out: "Eph", "Gp130 GST-II", LOC646627 and LOC100289019. The Marimuthu et al. 2011 block is a **gastric** cancer study.
- **Tested run.** The full pipeline was run end to end on 14 Nov 2026 (R 4.6.2, 2 cores, 7 GB RAM; about 25 minutes after downloads). Results are in `m1_results.zip`.
- **Memory.** Script 01 reads only the two needed columns from each TCGA file instead of using `GDCprepare()`, which needs more than 8 GB RAM.
- **msigdbr.** Version 10+ downloads gene sets from zenodo.org on first use. Without that access, install the self-contained version 7.5.1:
  `install.packages("https://cloud.r-project.org/src/contrib/Archive/msigdbr/msigdbr_7.5.1.tar.gz", repos = NULL, type = "source")`
- **TCGA MSI** is available for only about 45% of tumours (TCGA 2012 subtype table), so by the pre-specified >50%-missing rule it is left out of the TCGA models; GSE39582 keeps it.
- **AI assistance.** Code drafted with an AI assistant; state this in the Methods.
