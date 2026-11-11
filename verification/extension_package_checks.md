# Extension-package verification

Date: November 11, 2026.

The public package passed all 216 checks in `extensions_validation_checks.csv`: 111 baseline checks, 30 joint-covariance checks, 52 ordinal-extension checks, 18 enriched-imputation checks, three public-source-manifest checks, and two module-scope checks. All six suites completed without recorded warnings. The separate `precision_retry_checks.csv` records 39 passing checks. These tests used the copied public R files in native Windows R 4.5.1. They include public-package fixtures without internal specifications or review notes.

The thirteen extension, integration, and numerical-retry R files are byte-identical to the frozen analysis sources. Original public R files and historical result files remain unchanged. `extension_source_map.csv` records the source and public SHA-256 values for copied and derived release artifacts. The final release manifest covers every Git-visible file except itself.

The previously verified combined priority 1–2 CLI reproduced four principal CSVs byte-for-byte: joint covariance summaries, covariance geometry, ordinal point estimates, and ordinal model status. `priority12_cli_comparison.json` and `priority12_cli_stage_status.csv` retain that evidence. Packaging did not refit the study outcome models.

The final enriched-imputation summaries use 500 completed chains: five scenarios, two seeds, and 50 chains per seed, each with 50 iterations. All 1,000 required outcome fits passed after four bounded numerical retries. The 996 initially valid fits were retained. The four original failures involved only the Q6 primary score tolerance. The public retry-status table preserves their initial and final numerical diagnostics; the archived initial fits remain unchanged. The 20 scenario/seed/outcome status rows account for every required fit.

All 25 coefficient Monte Carlo checks passed. The largest Monte Carlo error was 0.00340812346329626, below 0.01. The largest absolute difference between seeds was 0.0067209396218475, below 0.02. Every final scenario summary is available. The covariance files describe variation between completed-data estimates and do not supply pooled sampling uncertainty.

All 25 final scenario summary CSVs match the authoritative final assembly byte-for-byte. The published covariance summary changes only its four local absolute `source_run` values to `outputs/public_cc_20261111_v3`. Every other parsed cell is unchanged. The four-row retry-status table selects numerical-status columns from the initial failures. The scenario-status table aggregates the complete final fit registry.

The donor-support helper produced 20 scenario/seed/item rows from the selected generation snapshot. Independent calculations matched within 8.9e-16. An in-memory input with one missing-cell update moved to another iteration was rejected after the per-chain/iteration coverage check was added. Valid output remained unchanged. The helper exports no recipient or donor identifiers.

The public Figure A1 CLI completed using the published aggregate inputs. Its regenerated 32-row plot-data CSV is byte-identical to the published figure data. The saved publication PNG used Matplotlib 3.10.5; the public command check used Matplotlib 3.10.8. Visual inspection found all labels and bars legible and within the specified axes. The original PNG remains unchanged.

The release boundary audit inspected the Git-visible file set and every published CSV schema and string field. It found no raw questionnaires, worker identifiers, support-row tables, donor-row tables, model objects, imputed datasets, internal mathematical documents, manuscript files, review records, or source books/articles. Full imputation trace and jurisdiction-distribution files remain excluded because they include singleton groups. Local absolute paths are absent from the published aggregate files.

The first local precision-test invocation completed its test function but had a trailing shell/R-expression syntax error. The corrected invocation used a fresh directory and exited successfully with all 39 checks passing. Both local directories remain excluded from Git. No study model or input changed during this correction.

The staged whitespace check passed with `core.whitespace=blank-at-eol,blank-at-eof,space-before-tab,cr-at-eol`. This retains ordinary whitespace checks while recognizing the preserved native Windows CRLF endings. The repository's existing `* -text` rule keeps source bytes and SHA-256 values unchanged.
