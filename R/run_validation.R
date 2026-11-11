# Unified synthetic validation. Sourcing defines the wrapper without executing.
# After independent review:
# Rscript run_validation.R --output-dir NEW_VERIFICATION_DIRECTORY
run_validation <- function(script_dir, output_dir) {
  if (dir.exists(output_dir) || file.exists(output_dir)) stop("Use a new verification output directory; prior results are preserved.")
  if (!dir.create(output_dir, recursive = TRUE)) stop("Could not create verification directory.")
  dependencies <- c("data_prep.R", "model_engine.R", "imputation.R", "test_data_prep.R", "test_synthetic.R", "test_imputation.R")
  for (file in dependencies) source(file.path(script_dir, file), local = environment())
  suites <- list(
    data_prep = function() run_data_prep_tests(file.path(script_dir, "data_prep.R"), output_dir),
    model_engine = function() run_synthetic_tests(file.path(script_dir, "model_engine.R"), output_dir),
    imputation = function() run_imputation_tests())
  all_results <- list(); suite_status <- list()
  for (name in names(suites)) {
    warnings <- character()
    value <- tryCatch(withCallingHandlers(suites[[name]](), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
    failed <- inherits(value, "error")
    if (failed) {
      rows <- data.frame(test = "suite_execution", passed = FALSE, detail = conditionMessage(value))
      # Tests that report their own individual failures write CSV before raising.
      partial_file <- switch(name, data_prep = "data_prep_tests.csv", model_engine = "synthetic_tests.csv", "")
      if (nzchar(partial_file) && file.exists(file.path(output_dir, partial_file))) {
        partial <- read.csv(file.path(output_dir, partial_file), stringsAsFactors = FALSE)
        rows <- rbind(partial[, c("test", "passed", "detail")], rows)
      }
    } else {
      rows <- value
      if ("check" %in% names(rows)) names(rows)[names(rows) == "check"] <- "test"
      if (!"detail" %in% names(rows)) rows$detail <- ""
      rows <- rows[, c("test", "passed", "detail")]
    }
    rows$suite <- name
    all_results[[name]] <- rows
    suite_status[[name]] <- data.frame(suite = name, passed = !failed && all(rows$passed),
      recorded_checks = nrow(rows), error = if (failed) conditionMessage(value) else "",
      warnings = paste(unique(warnings), collapse = " | "))
    write.csv(rows, file.path(output_dir, paste0("validation_", name, ".csv")), row.names = FALSE)
    write.csv(do.call(rbind, suite_status), file.path(output_dir, "validation_suite_status.csv"), row.names = FALSE)
    write.csv(do.call(rbind, all_results), file.path(output_dir, "validation_checks.csv"), row.names = FALSE)
  }
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "validation_sessionInfo.txt"))
  source_files <- file.path(script_dir, c(dependencies, "run_validation.R"))
  write.csv(data.frame(file = normalizePath(source_files, winslash = "/"),
    md5 = unname(tools::md5sum(source_files))), file.path(output_dir, "validation_source_md5.csv"), row.names = FALSE)
  status <- do.call(rbind, suite_status)
  saveRDS(list(status = status, checks = all_results, session = sessionInfo()), file.path(output_dir, "validation_results.rds"))
  if (!all(status$passed)) stop("Synthetic validation failed in: ", paste(status$suite[!status$passed], collapse = ", "))
  message("All synthetic validation suites passed. Results: ", normalizePath(output_dir, winslash = "/"))
  invisible(status)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L || args[1] != "--output-dir") stop("Usage: Rscript run_validation.R --output-dir NEW_DIRECTORY")
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script) != 1L) stop("Could not identify validation script path.")
  run_validation(dirname(normalizePath(sub("^--file=", "", script))), args[2])
}
