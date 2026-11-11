# Priority 1 and 2 extensions. Sourcing defines functions without execution.
# Execute only after the recorded mathematical and implementation reviews.
# Rscript R/run_priority12.R --data input/data.csv \
#   --source-run outputs/public_cc_20261111_v3 --output-dir outputs/NEW_RUN

priority12_expected_input_sha256 <-
  "2DFCE56408F4411C01426556C250DBD932790CFD24293B278E63F32A8323EDED"

parse_priority12_args <- function(args) {
  options <- list(data = NULL, source_run = NULL, output_dir = NULL,
                  specification = NULL, review_record = NULL)
  if (length(args) %% 2L) stop("Arguments require --name value pairs.")
  if (length(args)) for (i in seq.int(1L, length(args), by = 2L)) {
    key <- gsub("-", "_", sub("^--", "", args[i]), fixed = TRUE)
    if (!startsWith(args[i], "--") || !key %in% names(options))
      stop("Unknown argument: ", args[i])
    if (!is.null(options[[key]])) stop("Repeated argument: ", args[i])
    options[[key]] <- args[i + 1L]
  }
  for (key in c("data", "source_run", "output_dir"))
    if (is.null(options[[key]]) || !nzchar(options[[key]]))
      stop("Missing required argument: --", gsub("_", "-", key, fixed = TRUE))
  options
}

run_priority12 <- function(args = commandArgs(trailingOnly = TRUE), script_dir = NULL) {
  if (is.null(script_dir)) {
    script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(script) != 1L) stop("Supply script_dir for an interactive call.")
    script_dir <- dirname(normalizePath(sub("^--file=", "", script), mustWork = TRUE))
  }
  script_dir <- normalizePath(script_dir, winslash = "/", mustWork = TRUE)
  opts <- parse_priority12_args(args)
  input <- normalizePath(opts$data, winslash = "/", mustWork = TRUE)
  source_run <- normalizePath(opts$source_run, winslash = "/", mustWork = TRUE)
  if (!dir.exists(source_run)) stop("--source-run must be a directory.")
  if (file.exists(opts$output_dir) || dir.exists(opts$output_dir))
    stop("Use a new output directory. Historical runs are never overwritten.")
  if (!dir.create(opts$output_dir, recursive = TRUE)) stop("Could not create the output directory.")
  outdir <- normalizePath(opts$output_dir, winslash = "/", mustWork = TRUE)
  started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  modules <- c("data_prep.R", "model_engine.R", "run_analysis.R",
               "joint_extensions.R", "ordinal_extensions.R")
  code_files <- file.path(script_dir, c(modules, "run_priority12.R"))
  environment_for_run <- new.env(parent = environment())
  stage_rows <- list()
  final_state <- "setup_started"
  on.exit({
    writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
    saveRDS(list(started_utc = started,
                 ended_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                 state = final_state, options = opts, stages = stage_rows),
            file.path(outdir, "run_state.rds"))
  }, add = TRUE)
  writeLines(c(paste("Started UTC:", started),
               paste("Working directory:", normalizePath(getwd(), winslash = "/")),
               "Arguments (one value per line):", args),
             file.path(outdir, "execution_command.txt"))
  saveRDS(opts, file.path(outdir, "run_options.rds"))
  write.csv(data.frame(time = character(), item = character(), explanation = character()),
            file.path(outdir, "deviation_log.csv"), row.names = FALSE)
  missing_code <- code_files[!file.exists(code_files)]
  if (length(missing_code)) stop("Required source files are missing: ", paste(missing_code, collapse = "; "))
  for (file in modules) source(file.path(script_dir, file), local = environment_for_run)

  # Hash the exact old fits used by the extension, not newly fitted substitutes.
  fit_ids <- unlist(lapply(c(6, 4), function(outcome)
    c(paste0("q", outcome, "_main_full"), paste0("q", outcome, "_omit_site_", 1:8))),
    use.names = FALSE)
  fit_files <- file.path(source_run, "fits", paste0(fit_ids, ".rds"))
  write.csv(data.frame(model_id = fit_ids, file = fit_files, present = file.exists(fit_files)),
            file.path(outdir, "source_fit_inventory.csv"), row.names = FALSE)
  optional_files <- character()
  for (key in c("specification", "review_record")) if (!is.null(opts[[key]])) {
    optional_files <- c(optional_files, normalizePath(opts[[key]], winslash = "/", mustWork = TRUE))
  }
  # In the internal package record the current specification when available.
  # Public reproducibility does not require unpublished research documents.
  default_spec <- file.path(dirname(script_dir), "specification", "mathematical_specification", "main.tex")
  if (is.null(opts$specification) && file.exists(default_spec)) optional_files <- c(optional_files, default_spec)
  raw_copy <- file.path(dirname(dirname(script_dir)), "dataset", "data.csv")
  raw_files <- if (file.exists(raw_copy)) raw_copy else character()
  provenance_files <- unique(c(input, code_files, fit_files[file.exists(fit_files)],
                               optional_files, raw_files))
  environment_for_run$write_sha256_manifest(provenance_files, outdir)
  provenance <- read.csv(file.path(outdir, "provenance_sha256.csv"), stringsAsFactors = FALSE,
                         check.names = FALSE)
  canonical <- function(path) tolower(normalizePath(path, winslash = "/", mustWork = TRUE))
  input_hash <- provenance$sha256[match(canonical(input), vapply(provenance$file, canonical, character(1L)))]
  if (length(input_hash) != 1L || is.na(input_hash) ||
      toupper(input_hash) != priority12_expected_input_sha256) {
    final_state <- "input_hash_failed"
    writeLines(c(paste("Expected:", priority12_expected_input_sha256),
                 paste("Observed:", paste(input_hash, collapse = "; "))),
               file.path(outdir, "input_validation_failure.txt"))
    stop("Input SHA-256 differs from the reviewed public release. No extension was run.")
  }
  if (length(raw_files)) {
    raw_hash <- provenance$sha256[match(canonical(raw_copy), vapply(provenance$file, canonical, character(1L)))]
    if (length(raw_hash) != 1L || is.na(raw_hash) || toupper(raw_hash) != toupper(input_hash))
      stop("The original dataset and analysis input copy have different SHA-256 hashes.")
  }
  data <- tryCatch(environment_for_run$read_study_data(input), error = function(e) {
    final_state <<- "input_validation_failed"
    writeLines(conditionMessage(e), file.path(outdir, "input_validation_failure.txt"))
    stop(e)
  })
  write.csv(data.frame(file = input, sha256 = input_hash, rows = nrow(data),
                       columns = ncol(attr(data, "raw_input")), specification_version = "1.4"),
            file.path(outdir, "input_manifest.csv"), row.names = FALSE)

  stages <- list(
    joint_extensions = function() environment_for_run$run_joint_extension(
      source_run, file.path(outdir, "joint_extensions")),
    ordinal_extensions = function() environment_for_run$run_ordinal_extensions(
      data, source_run, file.path(outdir, "ordinal_extensions")))
  result <- list()
  for (stage in names(stages)) {
    stage_started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    warnings <- character()
    value <- tryCatch(withCallingHandlers(stages[[stage]](), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
    failed <- inherits(value, "error")
    result[[stage]] <- if (failed) list(error = conditionMessage(value)) else value
    stage_rows[[stage]] <- data.frame(stage = stage, started_utc = stage_started,
      ended_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      execution_completed = !failed,
      error = if (failed) conditionMessage(value) else "",
      warnings = paste(unique(warnings), collapse = " | "), stringsAsFactors = FALSE)
    write.csv(do.call(rbind, stage_rows), file.path(outdir, "extension_stage_status.csv"), row.names = FALSE)
    saveRDS(result, file.path(outdir, "extension_results.rds"))
  }
  # A completed stage can contain planned model failures or unsupported summaries.
  # Its own registry is authoritative for scientific availability.
  status <- do.call(rbind, stage_rows)
  final_state <- if (all(status$execution_completed)) "execution_completed" else "stage_execution_failed"
  versions <- c("MASS", "ordinal")
  write.csv(data.frame(package = c("R", versions),
    version = c(as.character(getRversion()), vapply(versions, function(package)
      tryCatch(as.character(packageVersion(package)), error = function(e) "unavailable"), character(1L)))),
    file.path(outdir, "software_versions.csv"), row.names = FALSE)
  if (!all(status$execution_completed)) stop("Extension execution failed in: ",
    paste(status$stage[!status$execution_completed], collapse = ", "),
    ". All attempted outputs are retained in ", outdir)
  message("Priority 1 and 2 outputs: ", outdir,
          ". Inspect module registries for model and support failures.")
  invisible(result)
}

if (sys.nframe() == 0L) run_priority12()
