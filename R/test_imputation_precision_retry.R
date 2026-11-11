# Synthetic tests for the separate bounded numerical retry.
# Sourcing defines tests only. Do not execute before independent review.
# No study input, saved study fit, or imputation chain is read by this suite.
run_precision_retry_tests <- function(script_dir, fixture_dir) {
  script_dir <- normalizePath(script_dir, winslash = "/", mustWork = TRUE)
  if (dir.exists(fixture_dir) &&
      length(list.files(fixture_dir, all.files = TRUE, no.. = TRUE)))
    stop("Precision retry fixture directory must be new or empty.")
  dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)
  fixture_dir <- normalizePath(fixture_dir, winslash = "/", mustWork = TRUE)
  for (file in c("data_prep.R", "model_engine.R", "imputation_extensions.R",
                 "imputation_precision_retry.R"))
    source(file.path(script_dir, file), local = environment())
  results <- list()
  write_results <- function() {
    if (length(results)) write.csv(do.call(rbind, results),
      file.path(fixture_dir, "precision_retry_tests.csv"), row.names = FALSE)
  }
  on.exit(write_results(), add = TRUE)
  check <- function(name, condition, detail = "") {
    pass <- isTRUE(condition)
    results[[length(results) + 1L]] <<- data.frame(test = name, passed = pass,
      detail = as.character(detail), stringsAsFactors = FALSE)
    if (!pass) stop("Precision retry test failed: ", name)
    invisible(TRUE)
  }
  expect_error <- function(name, expression) {
    message <- tryCatch({ force(expression); NULL },
      error = function(e) conditionMessage(e))
    check(name, !is.null(message), if (is.null(message)) "No error." else message)
  }

  # These independent artificial questionnaire values give a full-rank design
  # and moderate four-level ordinal probabilities in all eight sites.
  completed <- mie_with_rng(code = {
    set.seed(831907L)
    n <- 800L
    d <- data.frame(Observation = as.character(seq_len(n)),
      Site = factor(rep(1:8, each = n / 8), levels = as.character(1:8)),
      Language = factor(rep(c("English", "Spanish"), length.out = n),
        levels = c("English", "Spanish")))
    for (variable in mie_variables())
      d[[variable]] <- sample(if (variable == "D2") 0:1 else 1:5,
        n, replace = TRUE)
    scores <- derive_scores(d, outcome = 6L, management = "stored")
    eta <- c(0, -.15, .1, -.05, .12, -.12, .08, .03)[as.integer(d$Site)] +
      .15 * scores$L_c + .1 * scores$E_c + .12 * scores$R_c -
      .1 * scores$M_c + .2 * scores$B_c + .25 * scores$A_c
    cumulative <- vapply(c(-1.1, -.1, 1.1),
      function(threshold) plogis(threshold - eta), numeric(n))
    level <- 1L + rowSums(runif(n) > cumulative)
    d$Q6 <- c(1L, 3L, 4L, 5L)[level]
    d
  })
  original <- completed
  original$Q6[seq(100L, 800L, by = 100L)] <- NA_integer_
  original$Q11[c(9L, 209L)] <- NA_integer_
  original$Q12[c(10L, 310L)] <- NA_integer_
  state <- list(version = "1.4", status = "complete",
    scenario = mie_scenarios()[1L, , drop = FALSE], seed = mie_seeds()[1L],
    chain = 1L, iteration = 50L, data = completed, failure = NULL)
  original_before <- original
  state_before <- state
  check("synthetic original passes existing variable validation",
    isTRUE(mie_validate(original)))
  mipr_validate_completed(original, state)

  analysis <- mie_analysis_data(original, state$data, 6L)
  scored <- derive_scores(analysis, outcome = 6L, management = "stored")
  predictors <- c("L_c", "E_c", "R_c", "M_c", "B_c", "A_c")
  prep <- .prepare_model(scored, "Y", predictors, "Site")
  prep_before <- prep
  id <- sprintf("enriched_baseline_seed%d_chain1_iter50_Q6", state$seed)
  reference <- fit_ordinal(scored, "Y", predictors, "Site", id)
  check("independent synthetic reference fit is accepted",
    reference$valid && reference$inference_valid, reference$status)
  check("synthetic observed outcome mask excludes eight completed outcomes",
    nrow(prep$X) == 792L && identical(prep$rows,
      as.character(original$Observation[!is.na(original$Q6)])))

  # A fixed small perturbation creates the sole primary score defect while
  # retaining the independently computed original BFGS alternate.
  theta <- reference$theta
  theta["B_c"] <- theta["B_c"] + 3e-5
  component <- ordinal_components(theta, prep$X, prep$y, prep$K)
  check("synthetic primary exceeds the original score limit",
    max(abs(component$score)) / nrow(prep$X) >= 1e-6)
  attempts <- reference$diagnostics$attempts
  attempts$primary <- list(solver = "synthetic fixed perturbation",
    convergence = 0L, theta = theta, logLik = component$logLik,
    max_score_per_n = max(abs(component$score)) / nrow(prep$X),
    warnings = character(), error = NULL)
  attempts$primary_refinement <- NULL
  attempts$selected_primary <- "primary"
  failed <- .finish_fit(prep, theta, component, id, attempts,
    "primary_score_tolerance", character(), "ordinal", reference$alternate)
  record <- list(scenario = "baseline", seed = state$seed, chain = 1L,
    iteration = 50L, outcome = 6L, n = nrow(prep$X), valid = FALSE,
    reason = failed$status, values = c(beta_B = NA_real_, beta_A = NA_real_,
      d = NA_real_), fit = failed)
  record_before <- record
  saveRDS(list(original = original, state = state, reference = reference,
    failed_record = record), file.path(fixture_dir, "synthetic_fixture.rds"),
    version = 3L)

  polished <- mipr_polish(theta, prep)
  check("Newton polish recovers the perturbed synthetic fit",
    polished$valid && max(abs(polished$comp$score)) / nrow(prep$X) < 1e-8,
    polished$reason)
  check("Newton polish retains independent parameter agreement",
    max(abs(polished$theta - reference$alternate)) <= .01 &&
      max(abs(polished$theta[c("B_c", "A_c")] -
        reference$alternate[c("B_c", "A_c")])) <= 1e-4)
  trace <- polished$trace
  accepted_trace <- trace[trace$accepted, , drop = FALSE]
  check("Newton trace records initial point and bounded steps",
    is.data.frame(trace) && nrow(trace) >= 2L && nrow(trace) <= 1551L &&
      nrow(accepted_trace) <= 51L &&
      identical(as.integer(accepted_trace$iteration),
        seq.int(0L, nrow(accepted_trace) - 1L)) &&
      all(trace$iteration %in% 0:50) &&
      all(table(trace$iteration[-1L]) <= 31L) &&
      all(trace$halving[-1L] %in% 0:30) &&
      all(trace$alpha[-1L] == 2^(-trace$halving[-1L])))
  check("Newton accepted steps preserve likelihood within roundoff",
    all(diff(accepted_trace$logLik) >= -100 * .Machine$double.eps *
      (1 + abs(head(accepted_trace$logLik, -1L)))))
  check("Newton polish leaves prepared data unchanged",
    isTRUE(polished$data_unchanged) && identical(prep, prep_before))

  rank_bad <- prep
  rank_bad$X[, "A_c"] <- rank_bad$X[, "B_c"]
  check("Newton rejects a rank-deficient design",
    !isTRUE(mipr_polish(theta, rank_bad)$valid))
  threshold_bad <- theta
  threshold_index <- ncol(prep$X) + seq_len(prep$K - 1L)
  threshold_bad[threshold_index[2L]] <- threshold_bad[threshold_index[1L]] + 1e-8
  check("Newton rejects a threshold gap below the original limit",
    !isTRUE(mipr_polish(threshold_bad, prep)$valid))
  nonfinite <- theta
  nonfinite[1L] <- Inf
  check("Newton rejects a nonfinite initial parameter",
    !isTRUE(mipr_polish(nonfinite, prep)$valid))

  retry <- mipr_retry_record(record, original, state, expected_n = 792L)
  check("sole score failure is accepted after both fixed starts",
    retry$eligible && retry$accepted && isTRUE(retry$record$valid) &&
      length(retry$polishes) == 2L &&
      all(vapply(retry$polishes, function(x) isTRUE(x$valid), logical(1L))),
    retry$reason)
  check("retry preserves the complete initial failed record",
    identical(retry$initial_record, record_before) &&
      identical(record, record_before))
  check("retry preserves the independent saved BFGS alternate",
    identical(retry$record$fit$alternate, reference$alternate) &&
      identical(retry$record$fit$diagnostics$attempts$alternate,
        record$fit$diagnostics$attempts$alternate))
  check("retry retains exact observed rows, outcome, design, and eight sites",
    identical(retry$record$fit$rows, record$fit$rows) &&
      identical(retry$record$fit$outcome, record$fit$outcome) &&
      identical(retry$record$fit$design, record$fit$design) &&
      retry$record$fit$G == 8L && retry$record$fit$n == 792L)
  check("retry leaves original and completed questionnaire data unchanged",
    identical(original, original_before) && identical(state, state_before))

  other_failure <- record
  other_failure$fit$status <- "primary_score_tolerance;invalid_probabilities"
  other_failure$fit$diagnostics$reasons <- c("primary_score_tolerance", "invalid_probabilities")
  other_failure$reason <- other_failure$fit$status
  excluded <- mipr_retry_record(other_failure, original, state, expected_n = 792L)
  check("multiple failures are excluded from the retry",
    !isTRUE(excluded$eligible) && !isTRUE(excluded$accepted) &&
      identical(excluded$initial_record, other_failure))
  investigated <- record
  investigated$fit$diagnostics$investigation <- "coefficient_growth_or_alternate_instability"
  excluded <- mipr_retry_record(investigated, original, state, expected_n = 792L)
  check("an investigation requirement excludes the retry",
    !isTRUE(excluded$eligible) && !isTRUE(excluded$accepted))
  bad_alternate <- record
  bad_alternate$fit$diagnostics$attempts$alternate$convergence <- 1L
  excluded <- mipr_retry_record(bad_alternate, original, state, expected_n = 792L)
  check("an unaccepted original alternate cannot validate a retry",
    !isTRUE(excluded$accepted))

  wrong_data <- state
  wrong_data$data$Q11[1L] <- if (wrong_data$data$Q11[1L] == 5L) 4L else 5L
  expect_error("changed observed questionnaire cells are rejected",
    mipr_validate_completed(original, wrong_data))
  expect_error("a changed observed outcome mask is rejected",
    mipr_retry_record(record, original, state, expected_n = 793L))
  changed_mask <- original
  changed_mask$Q6[1L] <- NA_integer_
  changed_mask$Q6[100L] <- completed$Q6[100L]
  expect_error("an equal-sized changed observed outcome mask is rejected",
    mipr_retry_record(record, changed_mask, state, expected_n = 792L))
  wrong_design <- record
  wrong_design$fit$design[1L, "B_c"] <- wrong_design$fit$design[1L, "B_c"] + .1
  expect_error("saved design mismatch is rejected",
    mipr_retry_record(wrong_design, original, state, expected_n = 792L))
  wrong_id <- record
  wrong_id$fit$model_id <- "different_model"
  expect_error("saved model identifier mismatch is rejected",
    mipr_retry_record(wrong_id, original, state, expected_n = 792L))
  wrong_outcome <- record
  wrong_outcome$fit$outcome[1L] <- (wrong_outcome$fit$outcome[1L] %% 4L) + 1L
  expect_error("saved outcome mismatch is rejected",
    mipr_retry_record(wrong_outcome, original, state, expected_n = 792L))
  wrong_rows <- record
  wrong_rows$fit$rows <- rev(wrong_rows$fit$rows)
  expect_error("saved row order mismatch is rejected",
    mipr_retry_record(wrong_rows, original, state, expected_n = 792L))

  # These cloned identities test registry assembly only. They do not represent
  # additional completed chains or independently fitted synthetic outcomes.
  second <- record
  second$chain <- 2L
  second$fit$model_id <- sub("chain1_", "chain2_", second$fit$model_id, fixed = TRUE)
  second_retry <- retry
  second_retry$initial_record <- second
  second_retry$record$chain <- 2L
  second_retry$record$fit$model_id <- second$fit$model_id
  valid_source <- retry$record
  valid_source$chain <- 3L
  valid_source$fit$model_id <- sub("chain1_", "chain3_", valid_source$fit$model_id, fixed = TRUE)
  registry <- list(record, second, valid_source)
  retry$index <- 1L
  second_retry$index <- 2L
  assembled <- mipr_overlay_registry(registry, list(retry, second_retry))
  check("assembly applies every accepted retry and retains valid sources",
    length(assembled) == length(registry) &&
      identical(assembled[[1L]], retry$record) &&
      identical(assembled[[2L]], second_retry$record) &&
      identical(assembled[[3L]], valid_source))
  expect_error("assembly rejects a caller-selected subset of eligible retries",
    mipr_overlay_registry(registry, list(retry)))
  expect_error("assembly rejects duplicate retry indices",
    mipr_overlay_registry(registry, list(retry, retry, second_retry)))
  changed_source <- retry
  changed_source$initial_record$n <- changed_source$initial_record$n + 1L
  expect_error("assembly rejects a retry bound to a different source",
    mipr_overlay_registry(registry, list(changed_source, second_retry)))
  excluded <- mipr_retry_record(other_failure, original, state, expected_n = 792L)
  excluded$index <- 4L
  retained <- mipr_overlay_registry(c(registry, list(other_failure)),
    list(retry, second_retry, excluded))
  check("assembly retains a noneligible invalid source without replacement",
    identical(retained[[4L]], other_failure) && !isTRUE(retained[[4L]]$valid))

  # Existing summary rules must retain complete chain groups and both seeds.
  summary_registry <- list()
  for (seed in mie_seeds()) for (outcome in c(6L, 4L)) for (chain in 1:3) {
    values <- c(beta_B = .1 + .01 * chain,
      beta_A = .2 + .01 * chain + .001 * (seed - mie_seeds()[1L]))
    if (outcome == 6L) values <- c(values, d = unname(values["beta_A"] - values["beta_B"]))
    summary_registry[[length(summary_registry) + 1L]] <- list(
      scenario = "baseline", seed = seed, chain = as.integer(chain),
      iteration = 50L, outcome = outcome, n = 792L, valid = TRUE,
      reason = "", values = values)
  }
  complete <- mie_summarize_fits(summary_registry, 3L, "baseline", 50L)
  check("complete synthetic registry permits both outcome summaries",
    all(complete$availability$status == "available") &&
      setequal(unique(complete$estimates$outcome), c(6L, 4L)))
  incomplete <- summary_registry
  incomplete[[3L]]$valid <- FALSE
  incomplete[[3L]]$values[] <- NA_real_
  unavailable <- mie_summarize_fits(incomplete, 3L, "baseline", 50L)
  check("one failed chain blocks the entire outcome across both seeds",
    !any(unavailable$estimates$outcome == 6L) &&
      any(unavailable$availability$status == "unavailable") &&
      all(unavailable$estimates$outcome == 4L))
  omitted <- mie_summarize_fits(summary_registry[-3L], 3L, "baseline", 50L)
  check("omitting a failed chain cannot produce a subset average",
    !any(omitted$estimates$outcome == 6L))
  duplicated <- summary_registry
  duplicated[[3L]]$chain <- 2L
  duplicate_summary <- mie_summarize_fits(duplicated, 3L, "baseline", 50L)
  check("duplicated chain identities cannot replace a missing chain",
    !any(duplicate_summary$estimates$outcome == 6L))
  # Public reproduction requires executable sources, not unpublished notes.
  source(file.path(script_dir, "run_imputation_extensions.R"), local = environment())
  source(file.path(script_dir, "run_imputation_precision_retry.R"), local = environment())
  required_sources <- c("data_prep.R", "model_engine.R", "imputation_extensions.R",
    "run_imputation_extensions.R", "test_imputation_extensions.R",
    "imputation_precision_retry.R", "run_imputation_precision_retry.R",
    "test_imputation_precision_retry.R")
  public_dir <- file.path(fixture_dir, "public_package", "R")
  dir.create(public_dir, recursive = TRUE)
  stopifnot(all(file.copy(file.path(script_dir, required_sources), public_dir)))
  public_snapshot <- mipr_code_snapshot(public_dir)
  check("public package works without an internal note or specification",
    nrow(public_snapshot) == length(required_sources) &&
      setequal(basename(public_snapshot$file), required_sources))
  note_dir <- file.path(dirname(public_dir), "specification", "reviews")
  dir.create(note_dir, recursive = TRUE)
  writeLines("Synthetic internal implementation note.",
    file.path(note_dir, "precision_retry_implementation_note.md"))
  noted_snapshot <- mipr_code_snapshot(public_dir)
  check("an available internal note is included in the binding",
    nrow(noted_snapshot) == length(required_sources) + 1L &&
      "precision_retry_implementation_note.md" %in% basename(noted_snapshot$file))
  incomplete_dir <- file.path(fixture_dir, "incomplete_package", "R")
  dir.create(incomplete_dir, recursive = TRUE)
  stopifnot(all(file.copy(file.path(script_dir, required_sources[-1L]), incomplete_dir)))
  expect_error("missing executable source remains fatal",
    mipr_code_snapshot(incomplete_dir))
  saveRDS(list(retry = retry, polished = polished),
    file.path(fixture_dir, "synthetic_retry_results.rds"), version = 3L)
  do.call(rbind, results)
}
