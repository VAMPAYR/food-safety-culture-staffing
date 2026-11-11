# Reproducing the three analysis extensions

The extensions examine joint staffing uncertainty, selected ordinal-model restrictions, and missing-predictor assumptions. They were developed after inspection of the original results. They preserve the original questions and estimates. [METHODS.md](METHODS.md) defines their targets and reporting limits.

## Setup and original models

Run from the repository root in native Windows PowerShell. Install the `ordinal` R package; the recorded ordinal run used version 2025.12.29 with R 4.5.1 and MASS 7.3-65. Every run saves its actual software versions. The provenance commands require Windows PowerShell and `certutil`.

Use the commands in this document for the public checkout. Frozen source headers retain some original analysis-directory examples. The public imputation command must include `--data input/data.csv`; the figure command must use the public script and aggregate directories shown in the README.

```powershell
$env:LC_ALL = 'C'
$env:LANG = 'C'
python scripts/fetch_data.py
Rscript R/run_priority123_validation.R --output-dir verification/local_extensions
Rscript R/run_analysis.R --stage cc --data input/data.csv --run-name local_cc
```

The input must match SHA-256 `2DFCE56408F4411C01426556C250DBD932790CFD24293B278E63F32A8323EDED`. The original complete-case run generates the two full models and sixteen jurisdiction-deletion models needed by the first extension. Those row-level model files are not distributed in Git.

The validation runner includes the original data-preparation, model-engine, and imputation suites plus the three extension suites. It reports each suite and retains failures. The optional argument `--include-baseline false` runs only the extension suites when the unchanged baseline has already passed. Use a new validation directory for a new run.

Validate the separate numerical-retry module before using it:

```powershell
Rscript -e "source('R/test_imputation_precision_retry.R'); run_precision_retry_tests('R', 'verification/local_precision_retry')"
```

## Joint covariance and ordinal models

```powershell
Rscript R/run_priority12.R --data input/data.csv --source-run outputs/local_cc --output-dir outputs/local_priority12
```

Use a new output directory. The runner records source-model hashes, code and input hashes, software versions, stage status, and failures. It does not require private planning or review files.

| Directory | Contents |
| --- | --- |
| `joint_extensions/` | Joint covariance summaries, covariance directions, jurisdiction contributions, reconstruction checks, and source-fit provenance. |
| `ordinal_extensions/` | Partial proportional-odds and five-category staffing fits, acceptance status, support counts, and supported probability summaries. |

Read the model and probability status files even when the runner completes. Unavailable models and unsupported comparisons retain their reasons. A completed script does not imply every planned result was estimable. The joint F(2,7) probabilities remain uncalibrated references. The new ordinal fits and probability summaries are point-only sensitivities.

The published results retain all four accepted ordinal sensitivity fits. No probability comparison met the support rule. Each staffing comparison retained two jurisdictions; the rule requires at least three. The Q6 direct comparison also retained only 14 workers, below the 30-worker minimum. The empty probability tables and their status files record these unavailable displays.

## Enriched imputation

The five scenario names are `baseline`, `both_lower`, `both_higher`, `busy_lower`, and `absence_lower`. The default `--scenario all` dispatches them separately. Each scenario has its own manifest, checkpoints, diagnostics, review, and extension counts. Four parallel workers are supported; independent streams make the planned chain output independent of task scheduling.

Generate the imputed data first:

```powershell
Rscript R/run_imputation_extensions.R --stage generate --output outputs/local_priority3 --data input/data.csv --scenario all --m 50 --iterations 50 --workers 4
Rscript R/run_imputation_extensions.R --stage diagnostics --output outputs/local_priority3 --scenario all
```

Generation stops before outcome fitting. For each scenario, inspect chain status, traces of imputed means and variances, category proportions, distributions by jurisdiction, donor reuse, between-item associations, matching-regression diagnostics, and staffing donor-set support. Record specific observations and the reason for accepting or extending the chains.

After accepting one scenario's exact saved diagnostics, record its decision and fit its outcome models. For example:

```powershell
Rscript R/run_imputation_extensions.R --stage record-review --output outputs/local_priority3 --scenario baseline --decision accept --reviewer 'Reviewer name' --note 'Specific diagnostic findings and the reason for acceptance.'
Rscript R/run_imputation_extensions.R --stage fit --output outputs/local_priority3 --scenario baseline --workers 4
Rscript R/run_imputation_extensions.R --stage summary --output outputs/local_priority3 --scenario baseline
```

Repeat the review and fitting steps for the other scenarios. A batch `--scenario all` invocation is available after every selected scenario has a valid recorded decision. Each review binds the exact generation and diagnostic fingerprints. It cannot be overwritten to change its decision. A changed generation needs a new diagnostic review.

### Iteration extensions

For unresolved drift or mixing, record `--decision extend_iterations` for the affected scenario, then extend both roots to the next iteration count from their saved states:

```powershell
Rscript R/run_imputation_extensions.R --stage record-review --output outputs/local_priority3 --scenario baseline --decision extend_iterations --reviewer 'Reviewer name' --note 'Specific unresolved trace pattern requiring more iterations.'
Rscript R/run_imputation_extensions.R --stage generate --output outputs/local_priority3 --scenario baseline --m 50 --iterations 100 --workers 4
```

Review the new diagnostics before fitting. If needed, repeat the same sequence from 100 to 200 iterations. Unresolved drift at 200 makes that scenario unavailable. A numerical chain failure can be reported separately without accepting failed chains. Do not change the prior, donor count, or predictor set to obtain an accepted result.

### Chain-count extensions

After outcome fitting, each scenario's summary checks all reported coefficient components across its two seeds. If the Monte Carlo rule requires more chains, extend that scenario from 50 to 100 and then, if required, to 200 chains. Preserve its selected iteration count:

```powershell
Rscript R/run_imputation_extensions.R --stage generate --output outputs/local_priority3 --scenario baseline --m 100 --iterations 50 --workers 4
```

This example applies only when the scenario remains at 50 iterations and its previous fitted summary requires the next chain count. For a scenario already extended to 100 or 200 iterations, supply that iteration count. Review the added chains before fitting. Do not increase iterations and chain counts in the same generation step.

Only workers with the outcome originally observed enter that outcome's analysis. Coefficients are averaged on the log scale before exponentiation. A failed required fit makes its scenario/outcome summary unavailable. Completed data sets are not replaced, and a successful subset is not averaged. These summaries have no pooled confidence intervals or hypothesis tests.

### Bounded numerical retry

The separate retry addresses a saved ordinal fit that fails only the primary score tolerance. It preserves the completed data, observed-outcome rows, design, likelihood, and original acceptance limits. It applies bounded Newton steps from the saved primary and alternate estimates, requires the original agreement checks, and targets a maximum absolute score divided by sample size below `1e-8`. Every failed initial attempt remains available locally. Ineligible or unsuccessful retries retain their unavailable status.

After all five scenarios have completed their initial fitting and summaries, use separate new directories for the retry and final assembly:

```powershell
Rscript R/run_imputation_precision_retry.R --stage retry --input outputs/local_priority3 --output outputs/local_priority3_retry
Rscript R/run_imputation_precision_retry.R --stage assemble --input outputs/local_priority3_retry --output outputs/local_priority3_final
```

The assembly requires the full original registry and binds the unchanged generation, diagnostics, reviews, code, input, and retry evidence. It replaces only accepted numerical retries. It then calls the original summary function on every required fit. The final directory contains each scenario's `summaries_m50_i50/` files for the published run. Inspect `availability.csv`, `monte_carlo_checks.csv`, and the model registries before using `estimates.csv`. Neither a retry nor assembly changes the imputation chain count or iteration count.

The published run retained 50 chains per seed and 50 iterations for all five scenarios. All 500 chains completed. Of the 1,000 required outcome fits, 996 passed initially and four Q6 fits passed the bounded retry. All final summaries are available, and all 25 coefficient Monte Carlo checks passed. `results/extensions/imputation/scenario_seed_status.csv` accounts for the complete registry; `numerical_retry_status.csv` records the four recovered fits.

### Aggregate donor-support diagnostics

The public diagnostic files include original missing-staffing counts within each observed-outcome population, chain completion status, and terminal imputed staffing means by seed. The terminal means use only originally missing staffing responses. Their differences compare each scenario with the enriched baseline for the same seed, outcome, and item.

Create the published donor-support summary from a run in which all five scenarios have 50 chains and 50 iterations:

```powershell
python scripts/summarize_donor_support.py --input-root outputs/local_priority3 --output-file outputs/local_donor_support.csv --m 50 --iterations 50
```

This helper requires the specified counts for every scenario and preserves existing output files. If all scenarios share another allowed count, supply that count explicitly. The public 20-row table groups by scenario, seed, and staffing item. It reports the fraction of donor sets with one represented response category, the mean number of represented categories, the mean donor response before and after weighting, their mean difference, and the mean change in the probability of selecting agreement (response 4 or 5), in percentage points. Each missing-cell update receives equal weight. These conditional donor-selection summaries differ from the achieved terminal imputed means.

The helper assumes the recorded seeds and the analyzed dataset's three missing Q11 and two missing Q12 responses. It checks complete chain/iteration update counts. It does not determine whether a selected snapshot is final or diagnostically accepted. Verify that status separately, especially when a scenario has extended beyond a retained earlier snapshot. The single-category fraction ranges from zero to one; agreement shifts use percentage points.

## Files retained locally

Local `outputs/` contains row-level models, imputation checkpoints, donor records, diagnostic fingerprints, and execution manifests. It is excluded from Git. Published results contain aggregate estimates, support status, and selected diagnostics only. Original results remain in their existing directories; the extension results have separate directories.

Full trace and jurisdiction-distribution files remain local because some groups contain one imputed response. The donor-support helper removes all recipient and donor identifiers before aggregation. The public diagnostic tables contain no recipient-level records.

When copying aggregate results for publication, remove local absolute path fields or replace them with documented relative source labels. Do not publish worker identifiers, support-row files, donor-row records, fitted-model RDS files, or imputed data.

The published `joint_covariance_summary.csv` replaces its local absolute `source_run` value with `outputs/public_cc_20261111_v3`. All other parsed cells remain unchanged. `verification/extension_source_map.csv` records the original and public SHA-256 values and identifies every metadata-only change. The relative source labels identify the archived analysis run; they do not require that historical directory to be distributed.
