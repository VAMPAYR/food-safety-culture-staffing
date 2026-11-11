# Focused runner scope checks. No study data or model estimation is needed.
run_extension_integration_tests <- function(script_dir, fixture_dir) {
  script_dir <- normalizePath(script_dir, winslash = "/", mustWork = TRUE)
  if (!dir.exists(fixture_dir)) dir.create(fixture_dir, recursive = TRUE)
  fixture_dir <- normalizePath(fixture_dir, winslash = "/", mustWork = TRUE)
  isolated <- new.env(parent = baseenv())
  for (file in c("data_prep.R", "model_engine.R", "ordinal_extensions.R"))
    sys.source(file.path(script_dir, file), envir = isolated)
  message_from <- function(expr) tryCatch({ force(expr); "no error" }, error = conditionMessage)
  # The pre-existing directory is intentionally rejected after dependencies pass.
  local_message <- message_from(isolated$run_ordinal_extensions(NULL, "unused", fixture_dir))
  stopifnot(identical(local_message, "Ordinal output directory must be new."))
  incomplete <- new.env(parent = baseenv())
  sys.source(file.path(script_dir, "ordinal_extensions.R"), envir = incomplete)
  missing_message <- message_from(incomplete$run_ordinal_extensions(NULL, "unused", fixture_dir))
  stopifnot(identical(missing_message, "Source data_prep.R and model_engine.R first."))
  data.frame(test = c("Dependencies resolve in an isolated module environment",
                      "Missing dependencies remain rejected"), passed = TRUE,
             detail = c(local_message, missing_message), stringsAsFactors = FALSE)
}
