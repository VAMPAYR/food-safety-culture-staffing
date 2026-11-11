# Staffing coverage and perceived coworker food safety practices

This repository contains R code, aggregate results, and figures for a secondary analysis of the public EHS-Net food-worker survey described by Kramer et al. (2023). The analysis asks whether busy-period and absence coverage contribute information about perceived sick-worker exclusion and handwashing after adjustment for four culture scores and jurisdiction.

The publication manuscript and private research documents are maintained separately. Survey records are obtained directly from CDC and are excluded from Git. No additional CDC data or clarification is needed to execute this analysis.

## Contents

- `R/`: input checks, score derivation, ordinal and binary models, clustered inference, missing-data sensitivity, diagnostics, and tests.
- `results/complete_case/`: aggregate complete-case outputs. The main analyses include 539 Q6 and 541 Q4 observations.
- `results/imputation/mi_estimates.csv`: point estimates from both 50-chain imputation runs.
- `results/extensions/`: aggregate summaries of joint covariance, ordinal-model sensitivity, and enriched missing-data analyses.
- `figures/`: publication figures and the aggregate rows used to draw them.
- `METHODS.md`: coding, estimation, diagnostics, and reporting choices needed to understand the code.
- `EXTENSIONS.md`: commands, diagnostic review, and output definitions for the three post-results extensions.
- `scripts/fetch_data.py`: verified public-data download, with no changes to parsed cell values.
- `verification/`: package checks and source provenance.

## Input

Data title: *Development of an Empirically Derived Measure of Food Safety Culture in Restaurants*.

Provider: Centers for Disease Control and Prevention. Dataset identifier: `37nu-tuw8`.

Catalog: https://catalog.data.gov/dataset/development-of-an-empirically-derived-measure-of-food-safety-culture-in-restaurants

CSV: https://data.cdc.gov/api/views/37nu-tuw8/rows.csv?accessType=DOWNLOAD

The analyzed file contains 579 records and 40 columns. Its SHA-256 is `2DFCE56408F4411C01426556C250DBD932790CFD24293B278E63F32A8323EDED`. The live CDC download was checked on November 11, 2026: all parsed cells and row order matched the analyzed file. CSV quoting differed. The downloader normalizes only serialization and requires the exact analyzed checksum. A later content change causes it to stop.

```powershell
python scripts/fetch_data.py
```

## Software and execution

The original analysis used Windows, R 4.5.1, MASS 7.3-65, and parallel 4.5.1. The ordinal-model extension also used `ordinal` 2025.12.29. Install `ordinal` before running the extension checks. Each run records its actual package versions. The data downloader uses Python 3's standard library. Appendix Figure A1 uses Python 3 and Matplotlib; the saved figure used Matplotlib 3.10.5. Provenance helpers use native Windows PowerShell, .NET, or `certutil` for SHA-256.

Run from the repository root. Put `Rscript` on PATH or replace it with the full path to your Rscript executable. These commands use fresh output directories and refuse to overwrite historical runs.

```powershell
$env:LC_ALL = 'C'
$env:LANG = 'C'
Rscript R/run_validation.R --output-dir verification/local_validation
Rscript R/test_input_integration.R --input input/data.csv --output-dir verification/local_input
Rscript R/run_analysis.R --stage cc --run-name local_cc
Rscript R/run_analysis.R --stage mi --run-name local_mi
Rscript R/inspect_imputation.R --run-dir outputs/local_mi
```

The imputation stage stops after generating both seed runs. Inspect every chain's status, iterations, means and variances over iterations, distributions by jurisdiction, donor reuse, and correlations. Write the diagnostic findings and acceptance rationale in a local Markdown file. The acceptance script checks data preservation and structure; it does not decide whether the chains have mixed adequately.

After that diagnostic inspection, record acceptance and resume:

```powershell
Rscript R/record_mi_diagnostic_acceptance.R outputs/local_mi local_review.md "Reviewer name"
Rscript R/run_analysis.R --stage mi --resume-mi outputs/local_mi --mi-review outputs/local_mi/review_m50_iter50.rds --mi-m 50 --mi-iterations 50
```

Persistent drift requires extending all chains to 100 and then 200 iterations. Excessive Monte Carlo variation requires extending both seed runs to 100 and then 200 chains. See `METHODS.md` and the runner's header for extension decisions. Preserve checkpoints and supply both current counts explicitly. Do not pool a selected subset of successful fits. Generated files include the source hashes, software versions, full model objects, diagnostics, and numerical failures.

Generate the two figures from the published aggregate results:

```powershell
Rscript R/make_publication_figures.R
```

To plot a new run, replace the corresponding aggregate result files only in a separate working copy after verifying their source and scope. The script checks the reported sample sizes and expected categories.

## Post-results extensions

Three additions examine the joint staffing covariance, selected ordinal-model restrictions, and missing-predictor assumptions. They were developed after inspection of the original results. They preserve the original questions and analyses.

The joint-covariance extension reuses saved original fits. Reproduce `outputs/local_cc` first, then run the extension checks and the combined covariance/ordinal runner:

```powershell
Rscript R/run_priority123_validation.R --output-dir verification/local_extensions
Rscript R/run_priority12.R --data input/data.csv --source-run outputs/local_cc --output-dir outputs/local_priority12
```

The enriched imputation has five scenarios and a separate diagnostic-review gate before outcome fitting. Follow [EXTENSIONS.md](EXTENSIONS.md) for generation, scenario-specific review, fitting, and continuation commands. Its coefficient summaries are point estimates. They have no pooled sampling intervals or hypothesis tests.

Generate Appendix Figure A1 from the published aggregate covariance summaries:

```powershell
python scripts/make_priority123_figures.py --input-dir results/extensions/joint --output-dir outputs/local_figure_A1
```

The script checks all 32 jurisdiction contributions before drawing the PNG and PDF. Its supplied directory arguments make the public command independent of the original analysis directory.

Saved fits and imputation checkpoints are generated locally and remain excluded from Git. The published package contains aggregate summaries and computational checks. It does not require unpublished research documents.

[Extension-package verification](verification/extension_package_checks.md) records the 216 extension/baseline checks, 39 numerical-retry checks, final fit counts, figure reproduction, and release-boundary audit.

## Interpretation

The outcomes concern workers' perceptions of coworker practices. Q12 covers any employee absence. Management-item coding history remains unresolved; the management coefficient is not interpreted directionally. The public file has jurisdiction codes but no restaurant linkage. The eight-jurisdiction inference and jackknife intervals are approximate; their nominal coverage is not established for this exact design.

The two staffing associations are tested jointly, and the Q6 absence-minus-busy comparison directly uses their covariance. The Q4 outcome is excluded from its own employee predictor score. The missing-data analysis reports point estimates without pooled sampling intervals or tests. One Q4 binary diagnostic remains unavailable because of separation and ill-conditioned information. No favorable replacement was substituted.

The questions were developed after earlier exploratory exposure to staffing associations. They are not prospectively registered confirmatory hypotheses. No causal or out-of-sample prediction claim is made.

## Source

Kramer A, Hoover ER, Hedeen N, et al. Development of an empirically derived measure of food safety culture in restaurants. *Journal of Food Protection*. 2023;86:100043. https://doi.org/10.1016/j.jfp.2023.100043

The dataset is distributed by CDC under its own terms. This repository does not assign a new license to the CDC data or redistribute methodological books and articles. A software license has not yet been selected.
