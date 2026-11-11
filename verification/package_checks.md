# Public-package verification

Date: November 11, 2026.

This record describes the original release. The subsequent extension release is documented in [extension_package_checks.md](extension_package_checks.md).

The package passed 111 synthetic checks: 44 data-preparation, 56 model-engine, and 11 imputation checks. All six input/descriptive/SHA-256 integration checks passed. No suite warning was recorded. The figure script completed its aggregate-value and percentage checks. These tests used the published code and native Windows R 4.5.1; no study model was refitted during packaging.

The public-data downloader reproduced the analyzed SHA-256 exactly. Parsed cell values and row order matched the original analysis input. A first live response failed the checksum, and the downloader stopped without writing analyzed input. A later response passed. The checksum requirement remains unchanged and rejects changed content.

The core model engine, score preparation, and imputation estimation code retain their analysis-run logic. Publication changes are limited to the default input path, the methods file included in provenance, the diagnostic helper's reviewer-name argument, a diagnostic title, and figure input/output directories. Historical checkpoint schema strings remain unchanged for compatibility.

These checks establish computational agreement within the tested scope. They do not establish nominal uncertainty coverage with eight jurisdictions or validate the missingness assumptions.
