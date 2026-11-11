# Synthetic validation of the priority 1-3 extensions.
# Sourcing this file never executes tests or reads the study data.
# Rscript R/run_priority123_validation.R --output-dir verification/NEW_DIRECTORY

run_priority123_validation <- function(script_dir, output_dir, include_baseline = TRUE) {
  script_dir <- normalizePath(script_dir, winslash = "/", mustWork = TRUE)
  if (file.exists(output_dir) || dir.exists(output_dir))
    stop("Use a new validation directory; prior checks are preserved.")
  if (!dir.create(output_dir, recursive = TRUE)) stop("Could not create validation directory.")
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  validation_env <- new.env(parent = environment())
  modules <- c("data_prep.R", "model_engine.R", "run_analysis.R", "imputation.R",
    "joint_extensions.R", "ordinal_extensions.R", "imputation_extensions.R",
    "test_joint_extensions.R", "test_ordinal_extensions.R", "test_imputation_extensions.R",
    "test_extension_integration.R")
  if (include_baseline) modules <- c(modules, "run_validation.R",
    "test_data_prep.R", "test_synthetic.R", "test_imputation.R")
  source_files <- file.path(script_dir, unique(c(modules, "run_priority123_validation.R")))
  missing <- source_files[!file.exists(source_files)]
  if (length(missing)) stop("Required validation files are missing: ", paste(missing, collapse = "; "))
  for (file in modules) source(file.path(script_dir, file), local = validation_env)
  validation_env$write_sha256_manifest(source_files, output_dir)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "validation_sessionInfo_start.txt"))

  suites <- list()
  if (include_baseline) suites$baseline <- function() {
    baseline_dir <- file.path(output_dir, "baseline")
    validation_env$run_validation(script_dir, baseline_dir)
    read.csv(file.path(baseline_dir, "validation_checks.csv"), stringsAsFactors = FALSE)
  }
  suites$joint_extensions <- function() validation_env$run_joint_extension_tests(
    file.path(script_dir, "joint_extensions.R"), file.path(output_dir, "joint_extensions"))
  suites$ordinal_extensions <- function() validation_env$run_ordinal_extension_tests(
    script_dir, file.path(output_dir, "ordinal_extensions"))
  suites$imputation_extensions <- function() validation_env$run_imputation_extension_tests(script_dir)
  suites$imputation_manifest <- function() validation_env$run_imputation_manifest_tests(
    script_dir, file.path(output_dir, "public_manifest_fixture"))
  suites$extension_integration <- function() validation_env$run_extension_integration_tests(
    script_dir, file.path(output_dir, "module_scope_fixture"))
  all_rows <- list(); suite_rows <- list()
  for (suite in names(suites)) {
    warnings <- character()
    started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    value <- tryCatch(withCallingHandlers(suites[[suite]](), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
    failed <- inherits(value, "error")
    if (!failed) {
      if (!is.data.frame(value)) {
        failed <- TRUE
        value <- simpleError("Test suite did not return a data frame of recorded checks.")
      } else {
        if ("check" %in% names(value) && !"test" %in% names(value))
          names(value)[names(value) == "check"] <- "test"
        if (!all(c("test", "passed") %in% names(value)) || !nrow(value) ||
            !is.logical(value$passed) || anyNA(value$passed)) {
          failed <- TRUE
          value <- simpleError("Test suite returned empty or invalid check records.")
        }
      }
    }
    if (failed) {
      rows <- data.frame(test = "suite_execution", passed = FALSE, detail = conditionMessage(value))
      partial_name <- switch(suite, baseline = "validation_checks.csv",
        joint_extensions = "joint_extension_tests.csv",
        ordinal_extensions = "ordinal_extension_tests.csv", "")
      partial_file <- file.path(output_dir, suite, partial_name)
      if (nzchar(partial_name) && file.exists(partial_file)) {
        partial <- read.csv(partial_file, stringsAsFactors = FALSE)
        if ("check" %in% names(partial) && !"test" %in% names(partial))
          names(partial)[names(partial) == "check"] <- "test"
        if (all(c("test", "passed") %in% names(partial))) {
          if (!"detail" %in% names(partial)) partial$detail <- ""
          rows <- rbind(partial[, c("test", "passed", "detail")], rows)
        }
      }
    } else {
      if (!"detail" %in% names(value)) value$detail <- ""
      rows <- value[, c("test", "passed", "detail")]
    }
    rows$suite <- suite
    all_rows[[suite]] <- rows
    suite_rows[[suite]] <- data.frame(suite = suite,
      passed = !failed && all(rows$passed), recorded_checks = nrow(rows),
      started_utc = started, ended_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      error = if (failed) conditionMessage(value) else "",
      warnings = paste(unique(warnings), collapse = " | "), stringsAsFactors = FALSE)
    write.csv(rows, file.path(output_dir, paste0("validation_", suite, ".csv")), row.names = FALSE)
    write.csv(do.call(rbind, all_rows), file.path(output_dir, "validation_checks.csv"), row.names = FALSE)
    write.csv(do.call(rbind, suite_rows), file.path(output_dir, "validation_suite_status.csv"), row.names = FALSE)
  }
  status <- do.call(rbind, suite_rows)
  saveRDS(list(status = status, checks = all_rows, session = sessionInfo(),
               include_baseline = include_baseline), file.path(output_dir, "validation_results.rds"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "validation_sessionInfo.txt"))
  if (!all(status$passed)) stop("Validation failed in: ", paste(status$suite[!status$passed], collapse = ", "))
  message("All priority 1-3 validation suites passed: ", output_dir)
  invisible(status)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) %in% c(2L, 4L) || args[1] != "--output-dir")
    stop("Usage: Rscript run_priority123_validation.R --output-dir NEW_DIRECTORY [--include-baseline true|false]")
  include_baseline <- TRUE
  if (length(args) == 4L) {
    if (args[3] != "--include-baseline" || !args[4] %in% c("true", "false"))
      stop("--include-baseline must be true or false.")
    include_baseline <- args[4] == "true"
  }
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script) != 1L) stop("Could not identify validation script path.")
  run_priority123_validation(dirname(normalizePath(sub("^--file=", "", script), mustWork = TRUE)),
                            args[2], include_baseline)
}
