# Run only after the recorded implementation review and synthetic checks.
# Rscript run_analysis.R --stage cc|mi|all [--run-name NAME]
# Resume MI only: --stage mi --resume-mi EXISTING_RUN --mi-review REVIEW.rds
# Review file schema, created after inspection of saved chain diagnostics:
# saveRDS(list(reviewer="Reviewer name", notes="Diagnostic findings and rationale",
#              decision="accept"), "review.rds")
# Other decisions are extend_iterations, extend_chains, and unavailable.
# Extension examples append --mi-iterations 100 or --mi-m 100. Preserve the
# other current count explicitly if it differs from 50. Extensions stop for
# review of the newly generated terminal chains before fitting them.
# Sourcing defines functions; it never runs an analysis.

analysis_predictors <- c("L_c", "E_c", "R_c", "M_c", "B_c", "A_c")

new_fit_registry <- function(outdir) {
  e <- new.env(parent = emptyenv())
  e$fits <- list()
  e$outdir <- outdir
  dir.create(file.path(outdir, "fits"), showWarnings = FALSE)
  e
}

register_fit <- function(registry, data, outcome, predictors, model_id, kind = "ordinal") {
  original_id <- model_id
  repeat_index <- 1L
  while (model_id %in% names(registry$fits) || file.exists(file.path(registry$outdir, "fits", paste0(model_id, ".rds")))) {
    repeat_index <- repeat_index + 1L
    model_id <- paste0(original_id, "_attempt", repeat_index)
  }
  result <- tryCatch({
    fitter <- if (kind == "ordinal") fit_ordinal else fit_binary
    fitter(data, outcome, predictors, site = "Site", model_id = model_id)
  }, error = function(e) list(model_id = model_id, coef = numeric(), V = NULL, valid = FALSE,
    inference_valid = FALSE, status = paste("pipeline caught:", conditionMessage(e)), logLik = NA_real_,
    n = nrow(data), G = length(unique(data$Site)), rows = data$Observation,
    probabilities = NULL, diagnostics = list(error = conditionMessage(e))))
  result$model_id <- model_id
  result$analysis_metadata <- list(outcome = outcome, predictors = predictors, kind = kind,
                                    row_identifiers = data$Observation)
  registry$fits[[model_id]] <- result
  saveRDS(result, file.path(registry$outdir, "fits", paste0(model_id, ".rds")))
  result
}

save_registry <- function(registry) {
  saveRDS(registry$fits, file.path(registry$outdir, "model_registry.rds"))
  rows <- lapply(registry$fits, function(f) data.frame(model_id = f$model_id, n = f$n, sites = f$G,
    valid = isTRUE(f$valid), inference_valid = isTRUE(f$inference_valid),
    status = paste(f$status, collapse = "; "), log_likelihood = f$logLik,
    stringsAsFactors = FALSE))
  if (length(rows)) write.csv(do.call(rbind, rows), file.path(registry$outdir, "model_failure_registry.csv"), row.names = FALSE)
}

scalar_result <- function(fit, weights, term, sensitivity = FALSE, point_only = FALSE) {
  if (!isTRUE(fit$valid) || !all(names(weights) %in% names(fit$coef))) {
    ans <- data.frame(model_id = fit$model_id, term = term, estimate = NA_real_, se = NA_real_,
      lower = NA_real_, upper = NA_real_, exp_estimate = NA_real_, exp_lower = NA_real_, exp_upper = NA_real_,
      statistic = NA_real_, df1 = NA_real_, df2 = NA_real_, p_reference = NA_real_, valid = FALSE,
      status = paste(fit$status, collapse = "; "))
  } else {
    ans <- contrast_summary(fit, weights, label = term)
  }
  if (sensitivity) ans$p_reference <- NA_real_
  if (point_only) {
    for (nm in intersect(c("se", "lower", "upper", "exp_lower", "exp_upper", "statistic", "df1", "df2", "p_reference"), names(ans))) ans[[nm]] <- NA_real_
    # Point validity is independent of the intentionally excluded interval.
    ans$valid <- isTRUE(fit$valid)
    ans$status <- if (isTRUE(fit$valid)) "point estimate only; site deletion" else paste(fit$status, collapse = "; ")
  }
  ans
}

staffing_results <- function(fit, outcome, sensitivity = TRUE, point_only = FALSE,
                             terms = c("B_c", "A_c")) {
  rows <- lapply(terms, function(term) scalar_result(fit, setNames(1, term), term, sensitivity, point_only))
  if (outcome == 6 && identical(terms, c("B_c", "A_c"))) {
    rows[[length(rows) + 1L]] <- scalar_result(fit, c(A_c = 1, B_c = -1), "A_minus_B", sensitivity, point_only)
  }
  ans <- do.call(rbind, rows)
  ans$outcome <- outcome
  ans$n <- fit$n
  ans
}

add_change_from <- function(results, reference) {
  vals <- reference$estimate[match(results$term, reference$term)]
  results$reference_model <- reference$model_id[1]
  results$reference_estimate <- vals
  results$change_from_reference <- results$estimate - vals
  results
}

write_table <- function(rows, outdir, name) {
  if (is.data.frame(rows)) value <- rows else if (length(rows)) value <- do.call(rbind, rows) else value <- data.frame(status = "No rows produced; inspect model registry.")
  write.csv(value, file.path(outdir, paste0(name, ".csv")), row.names = FALSE, na = "")
}

write_sha256_manifest <- function(files, outdir) {
  # Literal PowerShell strings preserve spaces and punctuation in Windows paths.
  ps_literal <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
  files <- normalizePath(files, winslash = "\\", mustWork = TRUE)
  destination <- file.path(normalizePath(outdir, winslash = "\\"), "provenance_sha256.csv")
  script <- file.path(outdir, "provenance_hash.ps1")
  writeLines(c("$ErrorActionPreference = 'Stop'", "$sourceFiles = @(",
    paste0("  ", ps_literal(files), ifelse(seq_along(files) < length(files), ",", "")), ")",
    "$manifest = @('file,bytes,sha256')",
    "foreach ($sourceFile in $sourceFiles) {",
    "  $stream = $null; $sha = $null",
    "  try {",
    "    $stream = [System.IO.File]::OpenRead($sourceFile)",
    "    $sha = [System.Security.Cryptography.SHA256]::Create()",
    "    $bytes = $stream.Length",
    "    $digest = [System.BitConverter]::ToString($sha.ComputeHash($stream)).Replace('-', '')",
    "  } finally {",
    "    if ($null -ne $stream) { $stream.Dispose() }",
    "    if ($null -ne $sha) { $sha.Dispose() }",
    "  }",
    "  $escapedPath = $sourceFile.Replace('\"', '\"\"')",
    "  $manifest += ('\"{0}\",{1},\"{2}\"' -f $escapedPath, $bytes, $digest)",
    "}",
    paste0("[System.IO.File]::WriteAllLines(", ps_literal(destination),
           ", [string[]]$manifest, [System.Text.UTF8Encoding]::new($false))")), script)
  shell <- Sys.which("powershell.exe")
  if (!nzchar(shell)) shell <- Sys.which("pwsh.exe")
  if (!nzchar(shell)) stop("Native PowerShell is required for the SHA-256 manifest.")
  status <- system2(shell, c("-NoProfile", "-NonInteractive", "-File", shQuote(normalizePath(script))),
                    stdout = file.path(outdir, "provenance_hash_stdout.txt"),
                    stderr = file.path(outdir, "provenance_hash_stderr.txt"))
  if (!identical(status, 0L) || !file.exists(destination)) stop("SHA-256 manifest generation failed; inspect provenance hash logs.")
  invisible(destination)
}

design_diagnostics <- function(data, outcome, outdir) {
  vif <- do.call(rbind, lapply(analysis_predictors, function(term) {
    others <- setdiff(analysis_predictors, term)
    f <- reformulate(c("Site", others), response = term)
    model <- lm(f, data = data, na.action = na.fail)
    y <- data[[term]]
    tss <- sum((y - mean(y))^2)
    rss <- sum(residuals(model)^2)
    data.frame(outcome = outcome, term = term, n = nrow(data),
      vif = if (tss <= 0) NA_real_ else if (rss <= .Machine$double.eps * tss) Inf else tss / rss,
      status = if (tss <= 0) "constant predictor" else "descriptive design diagnostic")
  }))
  write_table(vif, outdir, paste0("q", outcome, "_vif"))
  joint <- as.data.frame(table(Q11 = factor(data$Q11, levels = 1:5), Q12 = factor(data$Q12, levels = 1:5)))
  joint$outcome <- outcome
  write_table(joint, outdir, paste0("q", outcome, "_staffing_joint_frequencies"))
  corr <- data.frame(outcome = outcome, n = nrow(data),
    pearson = cor(data$B, data$A, method = "pearson"), spearman = cor(data$B, data$A, method = "spearman"))
  write_table(corr, outdir, paste0("q", outcome, "_staffing_correlations"))
}

run_complete_cases <- function(data, outdir, registry) {
  primary <- list(); joints <- list(); apparent <- list(); sensitivities <- list()
  po_slopes <- list(); po_probs <- list(); po_sites <- list(); site_delete <- list()
  jackknife <- list()
  main <- list()
  for (o in c(6, 4)) {
    cc <- complete_case_rows(data, o)
    d <- derive_scores(data[cc, , drop = FALSE], o)
    prefix <- paste0("q", o)
    full <- register_fit(registry, d, "Y", analysis_predictors, paste0(prefix, "_main_full"))
    reduced <- register_fit(registry, d, "Y", setdiff(analysis_predictors, c("B_c", "A_c")), paste0(prefix, "_main_reduced"))
    main[[prefix]] <- full
    baseline <- staffing_results(full, o, sensitivity = FALSE)
    baseline$p_reference[baseline$term != "A_minus_B"] <- NA_real_
    primary[[prefix]] <- baseline
    joint_label <- paste0("RQ", if (o == 6) 1 else 3)
    joints[[prefix]] <- if (isTRUE(full$valid)) joint_summary(full, c("B_c", "A_c"), label = joint_label) else
      scalar_result(full, c(B_c = 1), joint_label)
    joints[[prefix]]$df1 <- 2L
    same_rows <- identical(full$analysis_metadata$row_identifiers, reduced$analysis_metadata$row_identifiers)
    pair_check <- if (isTRUE(full$valid) && isTRUE(reduced$valid)) same_row_comparison(full, reduced) else NULL
    pair_valid <- same_rows && !is.null(pair_check) && isTRUE(pair_check$valid)
    gain <- if (pair_valid) full$logLik - reduced$logLik else NA_real_
    apparent[[prefix]] <- data.frame(outcome = o, n = nrow(d), same_rows = same_rows,
      full_log_likelihood = full$logLik, reduced_log_likelihood = reduced$logLik,
      delta_log_likelihood = gain, decrease_apparent_mean_negative_log_score = gain / nrow(d),
      status = if (!same_rows) "invalid: different rows" else if (is.na(gain)) "unavailable: invalid fit" else if (gain < -1e-6 * (1 + abs(full$logLik))) "invalid: negative nested gain" else "apparent same-row fit only")
    tryCatch(design_diagnostics(d, o, outdir), error = function(e) {
      write_table(data.frame(outcome = o, status = conditionMessage(e)), outdir, paste0(prefix, "_design_diagnostic_failure"))
    })

    # Every adjustment and staffing slope varies across these binary thresholds.
    for (k in 1:3) {
      dk <- d
      dk$Z <- as.integer(as.integer(d$Y) > k)
      binary <- register_fit(registry, dk, "Z", analysis_predictors, paste0(prefix, "_threshold_", k), "binary")
      ss <- staffing_results(binary, o)
      ss <- add_change_from(ss, baseline)
      sensitivities[[binary$model_id]] <- ss
      for (term in analysis_predictors) {
        rr <- scalar_result(binary, setNames(1, term), term, sensitivity = TRUE)
        rr$outcome <- o; rr$threshold <- k
        rr$ordinal_estimate <- if (isTRUE(full$valid) && term %in% names(full$coef)) full$coef[[term]] else NA_real_
        rr$change_from_ordinal <- rr$estimate - rr$ordinal_estimate
        po_slopes[[paste(prefix, k, term)]] <- rr
      }
      ordinal_p <- if (isTRUE(full$valid)) rowSums(full$probabilities[, (k + 1):4, drop = FALSE]) else rep(NA_real_, nrow(d))
      binary_p <- if (isTRUE(binary$valid)) as.numeric(binary$probabilities) else rep(NA_real_, nrow(d))
      discrepancy <- abs(binary_p - ordinal_p)
      available <- isTRUE(full$valid) && isTRUE(binary$valid)
      po_probs[[paste(prefix, k)]] <- data.frame(outcome = o, threshold = k, n = nrow(d),
        mean_absolute_probability_difference = if (available) mean(discrepancy) else NA_real_,
        maximum_absolute_probability_difference = if (available) max(discrepancy) else NA_real_,
        status = if (available) "diagnostic" else "inconclusive: one or both fits invalid")
      for (s in 1:8) {
        idx <- as.character(d$Site) == as.character(s)
        po_sites[[paste(prefix, k, s)]] <- data.frame(outcome = o, threshold = k, site = s, n = sum(idx),
          observed_above = if (any(idx)) mean(dk$Z[idx]) else NA_real_,
          ordinal_fitted_above = if (any(idx) && isTRUE(full$valid)) mean(ordinal_p[idx]) else NA_real_,
          binary_fitted_above = if (any(idx) && isTRUE(binary$valid)) mean(binary_p[idx]) else NA_real_)
      }
    }

    dc <- d
    dc$B_category <- staffing_factor(dc$B)
    dc$A_category <- staffing_factor(dc$A)
    category <- register_fit(registry, dc, "Y", c("L_c", "E_c", "R_c", "M_c", "B_category", "A_category"), paste0(prefix, "_staffing_categories"))
    category_terms <- c("B_categorydisagreement", "B_categoryagreement", "A_categorydisagreement", "A_categoryagreement")
    cat_summary <- staffing_results(category, o, terms = category_terms)
    cat_summary$reference_model <- NA_character_; cat_summary$reference_estimate <- NA_real_; cat_summary$change_from_reference <- NA_real_
    sensitivities[[category$model_id]] <- cat_summary

    for (variant in c("alternative_outcome", "drop_Q13", "reverse_Q13", "reverse_Q15", "reverse_Q17")) {
      dv <- if (variant == "alternative_outcome") {
        derive_scores(data[cc, , drop = FALSE], o, outcome_coding = if (o == 6) "five" else "three")
      } else derive_scores(data[cc, , drop = FALSE], o, management = variant)
      vf <- register_fit(registry, dv, "Y", analysis_predictors, paste0(prefix, "_", variant))
      sensitivities[[vf$model_id]] <- add_change_from(staffing_results(vf, o), baseline)
    }

    cov_complete <- complete.cases(d[, c("D1", "D2", "D5", "Language"), drop = FALSE])
    dw <- d[cov_complete, , drop = FALSE]
    subset_fit <- register_fit(registry, dw, "Y", analysis_predictors, paste0(prefix, "_worker_subset_base"))
    adjusted <- register_fit(registry, dw, "Y", c(analysis_predictors, "D1", "D5", "D2", "Language"), paste0(prefix, "_worker_adjusted"))
    subset_summary <- staffing_results(subset_fit, o)
    sensitivities[[subset_fit$model_id]] <- add_change_from(subset_summary, baseline)
    sensitivities[[adjusted$model_id]] <- add_change_from(staffing_results(adjusted, o), subset_summary)

    deleted_fits <- setNames(vector("list", 8L), as.character(1:8))
    for (site in 1:8) {
      dl <- d[as.character(d$Site) != as.character(site), , drop = FALSE]
      dl$Site <- factor(as.character(dl$Site), levels = as.character(setdiff(1:8, site)))
      lf <- register_fit(registry, dl, "Y", analysis_predictors, paste0(prefix, "_omit_site_", site))
      deleted_fits[[as.character(site)]] <- lf
      lsum <- add_change_from(staffing_results(lf, o, point_only = TRUE), baseline)
      lsum$omitted_site <- site
      site_delete[[lf$model_id]] <- lsum
    }
    save_registry(registry)
    # The scalar CV3 sensitivity uses the full fit and all eight deletions.
    # The engine rejects the entire outcome sensitivity if any required fit is
    # invalid. No successful-fit subset or additional testing family is used.
    jackknife[[prefix]] <- cluster_jackknife_staffing(full, deleted_fits, paste0("Q", o))
  }
  primary_table <- do.call(rbind, primary)
  joint_table <- do.call(rbind, joints)
  secondary_raw <- c(RQ2 = primary_table$p_reference[primary_table$outcome == 6 & primary_table$term == "A_minus_B"],
                     RQ3 = joint_table$p_reference[joint_table$term == "RQ3"])
  if (length(secondary_raw) != 2L) stop("Secondary family must contain exactly RQ2 and RQ3.")
  secondary <- data.frame(question = c("RQ2", "RQ3"), raw_reference_p = unname(secondary_raw),
    holm_reference_p = if (all(is.finite(secondary_raw))) p.adjust(secondary_raw, method = "holm") else rep(NA_real_, 2),
    status = if (all(is.finite(secondary_raw))) "two-member secondary family; approximate reference only" else "family adjustment unavailable: at least one reference calculation failed")
  write_table(primary_table, outdir, "main_staffing_estimates")
  write_table(joint_table, outdir, "main_joint_reference_tests")
  write_table(secondary, outdir, "secondary_holm_reference")
  write_table(apparent, outdir, "apparent_fit_gains")
  write_table(sensitivities, outdir, "sensitivity_staffing_estimates")
  write_table(po_slopes, outdir, "proportional_odds_six_slopes")
  write_table(po_probs, outdir, "proportional_odds_probability_discrepancies")
  write_table(po_sites, outdir, "proportional_odds_site_observed_fitted")
  write_table(jackknife, outdir, "cluster_jackknife_staffing")
  saveRDS(jackknife, file.path(outdir, "cluster_jackknife_staffing.rds"))
  deletion <- do.call(rbind, site_delete)
  write_table(deletion, outdir, "site_deletion_point_estimates")
  ranges <- do.call(rbind, lapply(split(deletion, interaction(deletion$outcome, deletion$term, drop = TRUE)), function(z) {
    good <- is.finite(z$estimate)
    data.frame(outcome = z$outcome[1], term = z$term[1], attempted = nrow(z), valid = sum(good),
      minimum = if (any(good)) min(z$estimate[good]) else NA_real_, maximum = if (any(good)) max(z$estimate[good]) else NA_real_)
  }))
  write_table(ranges, outdir, "site_deletion_ranges")
  invisible(main)
}

parse_analysis_args <- function(args) {
  opts <- list(stage = "cc", run_name = NULL, resume_mi = NULL, mi_review = NULL,
               mi_m = 50L, mi_iterations = 50L, data = NULL)
  if (length(args) %% 2L != 0L) stop("Arguments require --name value pairs.")
  if (length(args)) for (i in seq(1L, length(args), by = 2L)) {
    key <- gsub("-", "_", sub("^--", "", args[i]), fixed = TRUE)
    if (!key %in% names(opts) || !startsWith(args[i], "--")) stop("Unknown argument: ", args[i])
    opts[[key]] <- args[i + 1L]
  }
  if (!opts$stage %in% c("cc", "mi", "all")) stop("--stage must be cc, mi, or all.")
  opts$mi_m <- as.integer(opts$mi_m); opts$mi_iterations <- as.integer(opts$mi_iterations)
  if (!opts$mi_m %in% c(50L, 100L, 200L) || !opts$mi_iterations %in% c(50L, 100L, 200L)) stop("MI chain and iteration counts must be 50, 100, or 200.")
  opts
}

run_analysis <- function(args = commandArgs(trailingOnly = TRUE), script_dir = NULL) {
  if (is.null(script_dir)) {
    script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(script) != 1L) stop("Supply script_dir when calling run_analysis interactively.")
    script_dir <- dirname(normalizePath(sub("^--file=", "", script)))
  }
  source(file.path(script_dir, "data_prep.R"), local = environment(run_complete_cases))
  source(file.path(script_dir, "model_engine.R"), local = environment(run_complete_cases))
  opts <- parse_analysis_args(args)
  package_dir <- dirname(script_dir)
  input <- if (is.null(opts$data)) file.path(package_dir, "input", "data.csv") else opts$data
  if (!is.null(opts$resume_mi)) {
    if (opts$stage != "mi") stop("--resume-mi requires --stage mi.")
    outdir <- normalizePath(opts$resume_mi, mustWork = TRUE)
    if (!file.exists(file.path(outdir, "mi_checkpoint.rds"))) stop("Resume directory has no MI checkpoint.")
    if (is.null(opts$mi_review)) stop("Resuming MI requires a recorded diagnostic review RDS file.")
  } else {
    name <- if (is.null(opts$run_name)) paste0("run_", format(Sys.time(), "%Y%m%dT%H%M%S")) else opts$run_name
    if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", name) || name %in% c(".", "..")) stop("Invalid run name.")
    outdir <- file.path(package_dir, "outputs", name)
    if (file.exists(outdir)) stop("Run directory already exists; historical runs are never overwritten.")
    if (!dir.create(outdir, recursive = TRUE)) stop("Could not create output directory.")
  }
  registry <- new_fit_registry(outdir)
  previous <- file.path(outdir, "model_registry.rds")
  if (!is.null(opts$resume_mi) && file.exists(previous)) registry$fits <- readRDS(previous)
  on.exit(save_registry(registry), add = TRUE)
  data <- tryCatch(read_study_data(input), error = function(e) {
    failure_name <- paste0("input_validation_failure_", format(Sys.time(), "%Y%m%dT%H%M%S"), ".txt")
    writeLines(conditionMessage(e), file.path(outdir, failure_name))
    stop(e)
  })
  if (is.null(opts$resume_mi)) {
    write_input_descriptives(data, outdir, input)
    writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
    write.csv(data.frame(package = c("R", "MASS", "parallel"),
      version = c(as.character(getRversion()), as.character(packageVersion("MASS")), as.character(packageVersion("parallel")))),
      file.path(outdir, "software_versions.csv"), row.names = FALSE)
    source_files <- c(input, list.files(script_dir, pattern = "\\.R$", full.names = TRUE),
                      file.path(package_dir, "METHODS.md"))
    source_files <- source_files[file.exists(source_files)]
    write.csv(data.frame(file = normalizePath(source_files, winslash = "/"), md5 = unname(tools::md5sum(source_files))),
      file.path(outdir, "provenance_md5.csv"), row.names = FALSE)
    write_sha256_manifest(source_files, outdir)
    write.csv(data.frame(time = character(), item = character(), explanation = character()), file.path(outdir, "deviation_log.csv"), row.names = FALSE)
    saveRDS(opts, file.path(outdir, "run_options.rds"))
  }
  if (opts$stage %in% c("cc", "all")) run_complete_cases(data, outdir, registry)
  if (opts$stage %in% c("mi", "all")) {
    source(file.path(script_dir, "imputation.R"), local = environment())
    reviewer <- if (is.null(opts$mi_review)) NULL else readRDS(opts$mi_review)
    callback <- function(data_completed, outcome, chain, seed) {
      d <- derive_scores(data_completed, outcome)
      id <- paste0("mi_seed", seed, "_chain", chain, "_iter", opts$mi_iterations, "_q", outcome)
      fit <- register_fit(registry, d, "Y", analysis_predictors, id)
      list(valid = isTRUE(fit$valid), coefficients = fit$coef[c("B_c", "A_c")], reason = fit$status, attempts = fit$diagnostics)
    }
    result <- run_mi(data, outdir, m = opts$mi_m, iterations = opts$mi_iterations,
                     fit_function = callback, diagnostic_review = reviewer)
    saveRDS(result, file.path(outdir, "mi_run_result.rds"))
  }
  save_registry(registry)
  message("Analysis output directory: ", normalizePath(outdir, winslash = "/"))
  invisible(outdir)
}

if (sys.nframe() == 0L) {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
  run_analysis(script_dir = script_dir)
}
