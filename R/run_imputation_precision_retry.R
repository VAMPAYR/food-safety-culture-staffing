# Native Windows driver. Definition-only on source; CLI explicitly selects retry
# or assemble. Run only after the independent implementation reviews clear.
mipr_source_environment <- function(script_dir, envir) {
  for (name in c("data_prep.R", "model_engine.R", "imputation_extensions.R",
                 "run_imputation_extensions.R", "imputation_precision_retry.R"))
    source(file.path(script_dir, name), local = envir)
  invisible(TRUE)
}

mipr_file_snapshot <- function(paths, role) {
  paths <- sort(unique(normalizePath(paths, winslash = "/", mustWork = TRUE)))
  hashes <- unname(tools::md5sum(paths))
  if (anyNA(hashes)) stop("Could not hash every provenance input.")
  data.frame(role = rep(role, length(paths)), file = paths, md5 = hashes,
    stringsAsFactors = FALSE)
}

mipr_run_snapshot <- function(run_dir) {
  mipr_file_snapshot(list.files(run_dir, recursive = TRUE, full.names = TRUE,
    all.files = TRUE, no.. = TRUE), "initial_run")
}

mipr_code_snapshot <- function(script_dir) {
  paths <- c(mie_source_paths(script_dir), file.path(script_dir,
    c("imputation_precision_retry.R", "run_imputation_precision_retry.R",
      "test_imputation_precision_retry.R")))
  # The public package omits internal review notes. All R sources remain
  # mandatory; local runs also bind the note when it is present.
  note <- file.path(dirname(script_dir), "specification", "reviews",
                    "precision_retry_implementation_note.md")
  if (file.exists(note)) paths <- c(paths, note)
  mipr_file_snapshot(paths, "reviewed_code_and_specification")
}

mipr_assert_snapshot <- function(snapshot) {
  if (any(!file.exists(snapshot$file)) ||
      !identical(unname(tools::md5sum(snapshot$file)), snapshot$md5))
    stop("A bound provenance input changed.")
  invisible(TRUE)
}

mipr_fresh_directory <- function(path, source_root) {
  if (file.exists(path) || dir.exists(path)) stop("Output path must be new.")
  parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
  resolved <- paste0(parent, "/", basename(path))
  source_root <- normalizePath(source_root, winslash = "/", mustWork = TRUE)
  if (identical(tolower(resolved), tolower(source_root)) ||
      startsWith(tolower(resolved), paste0(tolower(source_root), "/")))
    stop("Output must be outside the source run.")
  if (!dir.create(resolved, recursive = FALSE)) stop("Could not create fresh output.")
  normalizePath(resolved, winslash = "/", mustWork = TRUE)
}

mipr_load_run <- function(run_dir, script_dir) {
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  source_paths <- mie_source_paths(script_dir)
  source_hashes <- data.frame(file = normalizePath(source_paths, winslash = "/"),
    sha256 = vapply(source_paths, mie_sha256, character(1L)), row.names = NULL)
  manifests <- originals <- list(); records <- list(); all_tasks <- list()
  expected_scenarios <- mie_scenarios()$scenario
  directories <- list.dirs(file.path(run_dir, "scenarios"), recursive = FALSE, full.names = FALSE)
  if (!setequal(directories, expected_scenarios)) stop("The complete five-scenario source run is required.")
  for (scenario in expected_scenarios) {
    directory <- file.path(run_dir, "scenarios", scenario)
    manifest_path <- file.path(directory, "manifest.rds"); manifest <- readRDS(manifest_path)
    if (!identical(manifest$version, "1.4") || !identical(manifest$scenario, scenario) ||
        !manifest$m %in% c(50L, 100L, 200L) || !manifest$iterations %in% c(50L, 100L, 200L) ||
        !identical(manifest$seeds, mie_seeds()) || !identical(manifest$source_hashes, source_hashes))
      stop("Generation specification, code, seeds, or identity changed: ", scenario)
    original_path <- file.path(directory, "original_encoded.rds")
    if (!identical(unname(tools::md5sum(original_path)), manifest$original_md5) ||
        !identical(mie_sha256(manifest$input), manifest$input_sha256))
      stop("Original input no longer matches its generation manifest.")
    original <- readRDS(original_path); mie_validate(original)
    if (sum(!is.na(original$Q6)) != 577L) stop("The original Q6 observed-outcome mask must contain 577 rows.")
    if (length(originals) && !identical(original, originals[[1L]]))
      stop("Scenarios do not share the exact original data.")
    diagnostic_path <- file.path(directory, sprintf("diagnostic_index_m%d_i%d.rds", manifest$m, manifest$iterations))
    review_path <- file.path(directory, sprintf("review_m%d_i%d.rds", manifest$m, manifest$iterations))
    mie_require_review(manifest_path, readRDS(diagnostic_path), review_path)
    generation_path <- file.path(directory, sprintf("generation_registry_m%d_i%d.csv", manifest$m, manifest$iterations))
    if (!file.exists(generation_path)) stop("Missing generation registry.")
    tasks <- mie_fit_tasks(directory, manifest$m, manifest$iterations, scenario)
    model_index <- read.csv(file.path(directory, sprintf("model_md5_m%d_i%d.csv", manifest$m, manifest$iterations)), stringsAsFactors = FALSE)
    if (anyNA(model_index$md5) || !identical(unname(tools::md5sum(model_index$file)), model_index$md5))
      stop("Initial model index no longer matches saved fits.")
    for (i in seq_len(nrow(tasks))) {
      task <- as.list(tasks[i, , drop = FALSE]); record <- mie_load_saved_fit(task)
      if (is.null(record)) stop("Full initial registry required; a fit is missing.")
      actual_path <- if (file.exists(paste0(task$fit_path, ".failure.rds"))) paste0(task$fit_path, ".failure.rds") else task$fit_path
      if (!actual_path %in% model_index$file) stop("Initial fit absent from its recorded model index.")
      task$source_path <- actual_path; task$index <- length(records) + 1L
      records[[task$index]] <- record; all_tasks[[task$index]] <- task
    }
    manifests[[scenario]] <- manifest; originals[[scenario]] <- original
  }
  list(run_dir = run_dir, manifests = manifests, originals = originals,
    records = records, tasks = all_tasks)
}

mipr_record_metrics <- function(record) {
  f <- record$fit
  value <- function(x) if (length(x) == 1L && is.finite(x)) unname(x) else NA_real_
  data.frame(valid = isTRUE(record$valid), reason = paste(record$reason, collapse = ";"),
    max_score_per_n = value(f$diagnostics$max_score_per_n),
    logLik = value(f$logLik), n = value(record$n), G = value(f$G),
    information_condition = value(f$diagnostics$information$condition),
    beta_B = value(f$coef["B_c"]), beta_A = value(f$coef["A_c"]), stringsAsFactors = FALSE)
}

run_precision_retry <- function(run_dir, outdir, script_dir) {
  script_dir <- normalizePath(script_dir, winslash = "/", mustWork = TRUE)
  mipr_source_environment(script_dir, environment(run_precision_retry))
  source_snapshot <- mipr_run_snapshot(run_dir)
  loaded <- mipr_load_run(run_dir, script_dir)
  if (!identical(source_snapshot, mipr_run_snapshot(loaded$run_dir))) stop("Initial run changed while loading.")
  code_snapshot <- mipr_code_snapshot(script_dir)
  input_paths <- unique(vapply(loaded$manifests, function(x) x$input, character(1L)))
  input_snapshot <- mipr_file_snapshot(input_paths, "raw_input")
  outdir <- mipr_fresh_directory(outdir, loaded$run_dir)
  dir.create(file.path(outdir, "retries")); dir.create(file.path(outdir, "traces"))
  retries <- list(); retry_artifacts <- list(); flat <- list()
  for (index in seq_along(loaded$records)) {
    record <- loaded$records[[index]]; task <- loaded$tasks[[index]]
    retry_path <- retry_md5 <- NA_character_; result <- NULL
    if (!isTRUE(record$valid)) {
      state <- readRDS(task$path)
      result <- mipr_retry_record(record, loaded$originals[[record$scenario]], state,
        expected_n = sum(!is.na(loaded$originals[[record$scenario]][[paste0("Q", record$outcome)]])))
      result$index <- index
      result$provenance <- mipr_file_snapshot(c(task$source_path, task$path), "source_fit_and_completed_chain")
      retry_path <- file.path(outdir, "retries", basename(task$fit_path))
      saveRDS(result, retry_path); retry_md5 <- unname(tools::md5sum(retry_path))
      retries[[length(retries) + 1L]] <- result
      retry_artifacts[[length(retry_artifacts) + 1L]] <- data.frame(index = index,
        path = retry_path, md5 = retry_md5, accepted = result$accepted, stringsAsFactors = FALSE)
      if (length(result$polishes)) {
        traces <- lapply(names(result$polishes), function(start) {
          x <- result$polishes[[start]]$trace
          if (!nrow(x)) return(NULL)
          data.frame(start = start, x, stringsAsFactors = FALSE)
        })
        write.csv(mie_bind(Filter(Negate(is.null), traces)),
          file.path(outdir, "traces", sub("\\.rds$", ".csv", basename(task$fit_path))), row.names = FALSE)
      }
    }
    before <- mipr_record_metrics(record)
    after <- mipr_record_metrics(if (is.null(result)) record else result$record)
    names(before) <- paste0("initial_", names(before)); names(after) <- paste0("final_", names(after))
    flat[[index]] <- data.frame(index = index, scenario = record$scenario, seed = record$seed,
      chain = record$chain, iteration = record$iteration, outcome = record$outcome,
      source_path = task$source_path, source_md5 = unname(tools::md5sum(task$source_path)),
      chain_path = task$path, chain_md5 = unname(tools::md5sum(task$path)),
      retry_path = retry_path, retry_md5 = retry_md5,
      eligible = !is.null(result) && isTRUE(result$eligible),
      accepted_retry = !is.null(result) && isTRUE(result$accepted),
      retry_reason = if (is.null(result)) "initially_valid_retained" else result$reason,
      before, after, stringsAsFactors = FALSE)
  }
  # Completeness is checked before publishing any registry. No pooling occurs here.
  invisible(mipr_overlay_registry(loaded$records, retries))
  if (!identical(source_snapshot, mipr_run_snapshot(loaded$run_dir))) stop("Initial run changed during retry.")
  mipr_assert_snapshot(code_snapshot); mipr_assert_snapshot(input_snapshot)
  registry <- mie_bind(flat); registry_path <- file.path(outdir, "retry_registry.csv")
  write.csv(registry, registry_path, row.names = FALSE)
  provenance <- rbind(source_snapshot, code_snapshot, input_snapshot)
  write.csv(provenance, file.path(outdir, "provenance_md5.csv"), row.names = FALSE)
  trace_paths <- list.files(file.path(outdir, "traces"), full.names = TRUE)
  manifest <- list(version = "1", run_dir = loaded$run_dir, script_dir = script_dir,
    created = as.character(Sys.time()), protocol = mipr_protocol(),
    source_snapshot = source_snapshot, code_snapshot = code_snapshot, input_snapshot = input_snapshot,
    retries = mie_bind(retry_artifacts), registry_path = registry_path,
    registry_md5 = unname(tools::md5sum(registry_path)),
    trace_snapshot = mipr_file_snapshot(trace_paths, "retry_trace"),
    required_records = length(loaded$records), initial_invalid = sum(!registry$initial_valid),
    eligible = sum(registry$eligible), accepted = sum(registry$accepted_retry),
    remaining_invalid = sum(!registry$final_valid), summaries_written = FALSE)
  saveRDS(manifest, file.path(outdir, "retry_manifest.rds"))
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
  invisible(manifest)
}

assemble_precision_retry <- function(retry_dir, outdir, script_dir) {
  script_dir <- normalizePath(script_dir, winslash = "/", mustWork = TRUE)
  mipr_source_environment(script_dir, environment(assemble_precision_retry))
  retry_dir <- normalizePath(retry_dir, winslash = "/", mustWork = TRUE)
  manifest_path <- file.path(retry_dir, "retry_manifest.rds"); manifest <- readRDS(manifest_path)
  mipr_assert_snapshot(manifest$code_snapshot); mipr_assert_snapshot(manifest$input_snapshot)
  mipr_assert_snapshot(manifest$source_snapshot); mipr_assert_snapshot(manifest$trace_snapshot)
  if (!identical(manifest$code_snapshot, mipr_code_snapshot(script_dir)) ||
      !identical(manifest$source_snapshot, mipr_run_snapshot(manifest$run_dir)) ||
      !identical(unname(tools::md5sum(manifest$registry_path)), manifest$registry_md5) ||
      !identical(unname(tools::md5sum(manifest$retries$path)), manifest$retries$md5))
    stop("Retry code, registry, artifacts, or original run changed before assembly.")
  loaded <- mipr_load_run(manifest$run_dir, script_dir)
  if (length(loaded$records) != manifest$required_records) stop("Source registry count changed.")
  retries <- lapply(manifest$retries$path, readRDS)
  registry <- mipr_overlay_registry(loaded$records, retries)
  outdir <- mipr_fresh_directory(outdir, manifest$run_dir)
  if (startsWith(tolower(outdir), paste0(tolower(retry_dir), "/")))
    stop("Final summaries must be outside the retry directory.")
  results <- list()
  for (scenario in names(loaded$manifests)) {
    current <- loaded$manifests[[scenario]]
    directory <- file.path(outdir, "scenarios", scenario); dir.create(directory, recursive = TRUE)
    selected <- Filter(function(x) identical(x$scenario, scenario), registry)
    results[[scenario]] <- mie_write_fit_summary(selected, directory, current$m,
      current$iterations, scenarios = scenario)
  }
  provenance <- read.csv(manifest$registry_path, stringsAsFactors = FALSE)
  provenance$selected_path <- ifelse(provenance$accepted_retry, provenance$retry_path, provenance$source_path)
  provenance$selected_md5 <- ifelse(provenance$accepted_retry, provenance$retry_md5, provenance$source_md5)
  write.csv(provenance, file.path(outdir, "final_model_registry.csv"), row.names = FALSE)
  write.csv(rbind(manifest$source_snapshot, manifest$code_snapshot, manifest$input_snapshot,
    mipr_file_snapshot(c(manifest_path, manifest$registry_path, manifest$retries$path), "retry_evidence")),
    file.path(outdir, "provenance_md5.csv"), row.names = FALSE)
  saveRDS(list(retry_manifest_path = manifest_path, retry_manifest_md5 = unname(tools::md5sum(manifest_path)),
    required_records = length(registry), initial_invalid = manifest$initial_invalid,
    accepted_retries = manifest$accepted, remaining_invalid = sum(!vapply(registry, function(x) isTRUE(x$valid), logical(1L))),
    protocol = manifest$protocol, created = as.character(Sys.time())), file.path(outdir, "assembly_manifest.rds"))
  if (!identical(manifest$source_snapshot, mipr_run_snapshot(manifest$run_dir))) stop("Initial run changed during assembly.")
  mipr_assert_snapshot(manifest$code_snapshot); mipr_assert_snapshot(manifest$input_snapshot)
  invisible(results)
}

run_precision_retry_cli <- function(args = commandArgs(trailingOnly = TRUE), script_dir = NULL) {
  if (is.null(script_dir)) {
    script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(script) != 1L) stop("Supply script_dir.")
    script_dir <- dirname(normalizePath(sub("^--file=", "", script)))
  }
  options <- list(stage = NULL, input = NULL, output = NULL)
  if (length(args) %% 2L) stop("Use --name value pairs.")
  if (length(args)) for (i in seq(1L, length(args), 2L)) {
    name <- sub("^--", "", args[i])
    if (!startsWith(args[i], "--") || !name %in% names(options)) stop("Unknown option.")
    options[[name]] <- args[i + 1L]
  }
  if (any(vapply(options, is.null, logical(1L))) || !options$stage %in% c("retry", "assemble"))
    stop("Require --stage retry|assemble --input DIRECTORY --output NEW_DIRECTORY.")
  if (options$stage == "retry") run_precision_retry(options$input, options$output, script_dir)
  else assemble_precision_retry(options$input, options$output, script_dir)
}
if (sys.nframe() == 0L) run_precision_retry_cli()
