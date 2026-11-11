# Targeted public-input/descriptive/provenance integration check. No model fits.
run_input_integration <- function(script_dir, input, output_dir) {
  if (file.exists(output_dir)) stop("Use a new integration-test output directory.")
  if (!dir.create(output_dir, recursive = TRUE)) stop("Could not create integration-test output directory.")
  source(file.path(script_dir, "data_prep.R"), local = environment())
  source(file.path(script_dir, "run_analysis.R"), local = environment())
  warnings <- character()
  checks <- list()
  check <- function(name, pass) {
    checks[[length(checks) + 1L]] <<- data.frame(test = name, passed = isTRUE(pass))
  }
  error <- tryCatch(withCallingHandlers({
    d <- read_study_data(input)
    write_input_descriptives(d, output_dir, input)
    # A known SHA-256 vector verifies the .NET calculation independently of
    # external PowerShell modules. Write exactly three ASCII bytes, no newline.
    known_file <- file.path(output_dir, "sha256_known_vector_abc.txt")
    con <- file(known_file, open = "wb")
    tryCatch(writeBin(charToRaw("abc"), con), finally = close(con))
    files <- c(input, file.path(script_dir, c("data_prep.R", "run_analysis.R", "test_input_integration.R")), known_file)
    write_sha256_manifest(files, output_dir)
    manifest <- read.csv(file.path(output_dir, "provenance_sha256.csv"), colClasses = "character", check.names = FALSE)
    check("manifest_has_every_requested_file", nrow(manifest) == length(files) &&
      setequal(normalizePath(manifest$file, winslash = "/"), normalizePath(files, winslash = "/")))
    check("manifest_uses_full_sha256_digests", all(grepl("^[A-Fa-f0-9]{64}$", manifest$sha256)))
    known <- normalizePath(manifest$file, winslash = "/") == normalizePath(known_file, winslash = "/")
    check("known_sha256_vector_matches", sum(known) == 1L &&
      toupper(manifest$sha256[known]) == "BA7816BF8F01CFEA414140DE5DAE2223B00361A396177A9CB410FF61F20015AD")
    rq4 <- read.csv(file.path(output_dir, "rq4_six_cells.csv"))
    check("rq4_has_six_cells_and_four_levels", nrow(rq4) == 24L &&
      nrow(unique(rq4[, c("leadership_high", "absence_category")])) == 6L &&
      all(table(interaction(rq4$leadership_high, rq4$absence_category)) == 4L))
    groups <- split(rq4, interaction(rq4$leadership_high, rq4$absence_category))
    check("rq4_within_cell_outcome_counts_and_percentages", all(vapply(groups, function(z) {
      sum(z$outcome_n) == z$n[1] && if (z$n[1] > 0) abs(sum(z$outcome_percent) - 100) < 1e-8 else all(is.na(z$outcome_percent))
    }, logical(1))))
    check("no_input_descriptive_or_hash_warnings", length(warnings) == 0L)
    NULL
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) conditionMessage(e))
  if (!is.null(error)) check("integration_execution", FALSE)
  write.csv(do.call(rbind, checks), file.path(output_dir, "input_integration_checks.csv"), row.names = FALSE)
  writeLines(c(warnings, if (!is.null(error)) paste("ERROR:", error)), file.path(output_dir, "input_integration_messages.txt"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "input_integration_sessionInfo.txt"))
  if (!is.null(error)) stop(error)
  if (!all(vapply(checks, function(z) z$passed, logical(1)))) stop("Input integration checks failed.")
  invisible(do.call(rbind, checks))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 4L || args[1] != "--input" || args[3] != "--output-dir")
    stop("Usage: Rscript test_input_integration.R --input DATA.csv --output-dir NEW_DIRECTORY")
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  run_input_integration(dirname(normalizePath(sub("^--file=", "", script))), args[2], args[4])
}
