# Definition-only numerical retry for an unchanged completed-data ordinal fit.
# Existing model, imputation, acceptance, and summary implementations are reused.
mipr_protocol <- function() list(version = "1", max_steps = 50L,
  max_halvings = 30L, target_score_per_n = 1e-8,
  original_score_per_n = 1e-6, minimum_threshold_gap = 1e-6,
  maximum_parameter = 10, maximum_information_condition = 1e10,
  likelihood_roundoff_multiplier = 100)

mipr_is_score_failure <- function(record) {
  f <- record$fit
  !isTRUE(record$valid) && !is.null(f) && identical(f$valid, FALSE) &&
    identical(f$type, "ordinal") && identical(f$status, "primary_score_tolerance") &&
    identical(f$diagnostics$reasons, "primary_score_tolerance") &&
    !length(f$diagnostics$investigation)
}

mipr_validate_completed <- function(original, state) {
  if (!identical(state$status, "complete")) stop("Completed chain required.")
  completed <- state$data
  if (!is.data.frame(original) || !is.data.frame(completed) ||
      !identical(dim(original), dim(completed)) ||
      !identical(names(original), names(completed)) ||
      !identical(rownames(original), rownames(completed)) ||
      !identical(original$Observation, completed$Observation))
    stop("Completed-data identity or row order changed.")
  mie_validate(original); mie_validate(completed)
  for (nm in names(original)) {
    observed <- !is.na(original[[nm]])
    if (!identical(original[[nm]][observed], completed[[nm]][observed]))
      stop("Observed data changed: ", nm)
    if (!nm %in% mie_variables() && !identical(original[[nm]], completed[[nm]]))
      stop("Non-imputed data changed: ", nm)
  }
  if (anyNA(completed[, mie_variables(), drop = FALSE]))
    stop("Completed chain still contains missing imputation variables.")
  invisible(TRUE)
}

mipr_candidate <- function(theta, prep, score_limit = 1e-6) {
  protocol <- mipr_protocol(); p <- ncol(prep$X)
  bad <- character()
  if (length(theta) != p + prep$K - 1L || any(!is.finite(theta)))
    return(list(valid = FALSE, reasons = "nonfinite_or_mismatched_parameters", comp = NULL))
  if (max(abs(theta)) > protocol$maximum_parameter) bad <- c(bad, "absolute_parameter_above_10")
  gap <- min(diff(theta[p + seq_len(prep$K - 1L)]))
  if (!is.finite(gap) || gap < protocol$minimum_threshold_gap)
    return(list(valid = FALSE, reasons = c(bad, "threshold_gap_below_1e-6"), comp = NULL))
  comp <- ordinal_components(theta, prep$X, prep$y, prep$K)
  if (!is.finite(comp$logLik) || is.null(comp$score) || any(!is.finite(comp$score)))
    return(list(valid = FALSE, reasons = c(bad, "nonfinite_likelihood_or_score"), comp = comp))
  score <- max(abs(comp$score)) / nrow(prep$X)
  if (score >= score_limit) bad <- c(bad, "score_tolerance")
  probability_error <- max(abs(rowSums(comp$probabilities) - 1))
  if (any(!is.finite(comp$probabilities)) ||
      any(comp$probabilities < 0 | comp$probabilities > 1) ||
      !is.finite(probability_error) || probability_error > 1e-8)
    bad <- c(bad, "invalid_probabilities")
  info <- .matrix_check(comp$H)
  if (!info$valid || info$rank != length(theta)) bad <- c(bad, "information_not_positive_definite")
  if (info$condition > protocol$maximum_information_condition)
    bad <- c(bad, "information_condition_above_1e10")
  sandwich <- site_sandwich(comp$scores, comp$H, prep$site)
  if (!isTRUE(sandwich$valid)) bad <- c(bad, "site_covariance_unavailable")
  list(valid = !length(bad), reasons = unique(bad), comp = comp,
    metrics = list(max_score_per_n = score, logLik = comp$logLik,
      threshold_gap = gap, maximum_parameter = max(abs(theta)),
      information_condition = info$condition, information_min_eigen = info$min_eigen,
      information_rank = info$rank, probability_sum_error = probability_error),
    information = info)
}

mipr_agreement <- function(theta, other, prep) {
  first <- mipr_candidate(theta, prep); second <- mipr_candidate(other, prep)
  if (!first$valid || !second$valid)
    return(list(valid = FALSE, reasons = c(first$reasons, second$reasons)))
  staffing <- grep("^(B_c|A_c)", names(theta), value = TRUE)
  staff_difference <- if (length(staffing)) max(abs(theta[staffing] - other[staffing])) else 0
  parameter_difference <- max(abs(theta - other))
  likelihood_difference <- abs(first$comp$logLik - second$comp$logLik)
  se_ratio <- sqrt(diag(chol2inv(chol(second$comp$H))) /
    diag(chol2inv(chol(first$comp$H))))
  reasons <- character()
  if (staff_difference > 1e-4) reasons <- c(reasons, "staffing_alternate_disagreement")
  if (likelihood_difference > 1e-6 * (1 + abs(first$comp$logLik)))
    reasons <- c(reasons, "likelihood_alternate_disagreement")
  if (parameter_difference > 0.01) reasons <- c(reasons, "coefficient_growth_or_alternate_instability")
  if (any(!is.finite(se_ratio)) || max(se_ratio) > 1.1)
    reasons <- c(reasons, "standard_error_growth_under_tighter_tolerance")
  list(valid = !length(reasons), reasons = reasons,
    staff_difference = staff_difference, parameter_difference = parameter_difference,
    likelihood_difference = likelihood_difference, information_se_ratio = se_ratio)
}

mipr_polish <- function(theta, prep) {
  protocol <- mipr_protocol(); original_prep <- prep
  trace <- list(); iterates <- list(theta)
  result <- function(valid, reason, current) list(valid = valid, reason = reason,
    theta = theta, comp = current$comp, metrics = current$metrics,
    trace = mie_bind(trace), iterates = iterates,
    data_unchanged = identical(original_prep, prep), protocol = protocol)
  if (!is.matrix(prep$X) || any(!is.finite(prep$X)) ||
      qr(prep$X)$rank != ncol(prep$X))
    return(result(FALSE, "invalid_or_rank_deficient_design", list(comp = NULL)))
  current <- mipr_candidate(theta, prep, score_limit = Inf)
  if (!current$valid) return(result(FALSE, paste(current$reasons, collapse = ";"), current))
  trace[[1L]] <- data.frame(iteration = 0L, halving = NA_integer_, alpha = 0,
    logLik = current$comp$logLik, max_score_per_n = current$metrics$max_score_per_n,
    allowed_logLik_decrease = 0, accepted = TRUE, reason = "initial")
  for (iteration in seq_len(protocol$max_steps)) {
    if (current$metrics$max_score_per_n < protocol$target_score_per_n)
      return(result(TRUE, "target_reached", current))
    R <- chol(current$comp$H)
    direction <- drop(backsolve(R, forwardsolve(t(R), current$comp$score)))
    if (any(!is.finite(direction))) return(result(FALSE, "nonfinite_newton_direction", current))
    accepted <- FALSE
    allowance <- protocol$likelihood_roundoff_multiplier * .Machine$double.eps *
      (1 + abs(current$comp$logLik))
    for (halving in 0:protocol$max_halvings) {
      alpha <- 2^(-halving); candidate_theta <- theta + alpha * direction
      candidate <- mipr_candidate(candidate_theta, prep, score_limit = Inf)
      ll <- if (is.null(candidate$comp)) NA_real_ else candidate$comp$logLik
      accepted <- isTRUE(candidate$valid) && is.finite(ll) &&
        ll >= current$comp$logLik - allowance
      trace[[length(trace) + 1L]] <- data.frame(iteration = iteration,
        halving = halving, alpha = alpha, logLik = ll,
        max_score_per_n = if (is.null(candidate$metrics)) NA_real_ else candidate$metrics$max_score_per_n,
        allowed_logLik_decrease = allowance, accepted = accepted,
        reason = if (accepted) "accepted" else paste(c(candidate$reasons,
          if (isTRUE(candidate$valid)) "likelihood_decrease"), collapse = ";"))
      if (accepted) {
        theta <- candidate_theta; current <- candidate
        iterates[[length(iterates) + 1L]] <- theta
        break
      }
    }
    if (!accepted) return(result(FALSE, "line_search_exhausted", current))
  }
  result(current$metrics$max_score_per_n < protocol$target_score_per_n,
    if (current$metrics$max_score_per_n < protocol$target_score_per_n) "target_reached" else "step_cap_reached", current)
}

mipr_retry_record <- function(record, original, state, expected_n = NULL) {
  initial_record <- record
  rejected <- function(reason, eligible = FALSE, polishes = NULL, checks = NULL)
    list(eligible = eligible, accepted = FALSE, reason = reason,
      initial_record = initial_record, record = initial_record,
      polishes = polishes, checks = checks, protocol = mipr_protocol())
  mipr_validate_completed(original, state)
  if (!identical(record$scenario, state$scenario$scenario) || record$seed != state$seed ||
      record$chain != state$chain || record$iteration != state$iteration ||
      !record$outcome %in% c(4L, 6L)) stop("Saved fit and completed-chain identity mismatch.")
  if (!mipr_is_score_failure(record)) return(rejected("not_sole_primary_score_failure"))
  expected_id <- sprintf("enriched_%s_seed%d_chain%d_iter%d_Q%d", record$scenario,
    record$seed, record$chain, record$iteration, record$outcome)
  if (!identical(record$fit$model_id, expected_id)) stop("Saved model identifier changed.")
  completed <- mie_analysis_data(original, state$data, record$outcome)
  prepared <- derive_scores(completed, outcome = record$outcome, management = "stored")
  prep <- .prepare_model(prepared, "Y", c("L_c", "E_c", "R_c", "M_c", "B_c", "A_c"), "Site")
  f <- record$fit
  if (!identical(prep$X, f$design) || !identical(prep$y, f$outcome) ||
      !identical(prep$rows, f$rows) || !identical(prep$site, f$site) ||
      nrow(prep$X) != f$n || record$n != f$n ||
      (!is.null(expected_n) && f$n != expected_n) || f$G != 8L ||
      !identical(sort(as.character(prep$sites)), as.character(1:8)))
    stop("Saved design, outcome, observed-outcome rows, sample size, or sites changed.")
  if (!identical(names(f$theta), f$parameter_names) ||
      !identical(f$theta, f$coef) || !identical(f$theta, f$coefficients) ||
      !identical(names(f$alternate), names(f$theta)))
    return(rejected("saved_parameter_identity_mismatch"))
  selected <- f$diagnostics$attempts$selected_primary
  if (is.null(selected)) selected <- "primary"
  if (!isTRUE(f$diagnostics$attempts[[selected]]$convergence == 0) ||
      !identical(f$diagnostics$attempts[[selected]]$theta, f$theta) ||
      !isTRUE(f$diagnostics$attempts$alternate$convergence == 0) ||
      !identical(f$alternate, f$diagnostics$attempts$alternate$theta))
    return(rejected("saved_solver_convergence_or_alternate_identity_failed"))
  before <- mipr_candidate(f$theta, prep)
  alternate <- mipr_candidate(f$alternate, prep)
  if (!identical(before$reasons, "score_tolerance") || !alternate$valid)
    return(rejected("recomputed_source_has_other_numerical_defect", checks = list(before = before, alternate = alternate)))
  if (!identical(before$comp$H, f$H) || !identical(before$comp$scores, f$scores) ||
      !identical(before$comp$logLik, f$logLik))
    stop("Saved primary derivatives or likelihood differ on the reconstructed data.")
  # Agreement is assessed under the original limits even though primary score fails.
  source_comparison <- mipr_agreement_without_score(f$theta, f$alternate, prep)
  if (!source_comparison$valid) return(rejected("source_alternate_disagreement", checks = source_comparison))
  polishes <- list(primary = mipr_polish(f$theta, prep), alternate = mipr_polish(f$alternate, prep))
  if (!all(vapply(polishes, function(x) isTRUE(x$valid) && isTRUE(x$data_unchanged), logical(1L))))
    return(rejected("bounded_polish_not_accepted", TRUE, polishes))
  candidate <- polishes$primary$theta
  final <- mipr_candidate(candidate, prep, score_limit = mipr_protocol()$target_score_per_n)
  saved_alternate_check <- mipr_agreement(candidate, f$alternate, prep)
  second_start_check <- mipr_agreement(candidate, polishes$alternate$theta, prep)
  checks <- list(before = before$metrics, saved_alternate = alternate$metrics,
    after = final$metrics, original_agreement = source_comparison,
    saved_alternate_agreement = saved_alternate_check, second_start_agreement = second_start_check)
  if (!final$valid || !saved_alternate_check$valid || !second_start_check$valid)
    return(rejected("original_acceptance_gates_failed", TRUE, polishes, checks))
  attempts <- f$diagnostics$attempts
  attempts$precision_retry <- list(solver = "bounded natural-parameter Newton score polish",
    selected_start = "saved_primary", protocol = mipr_protocol(), polishes = polishes,
    checks = checks, convergence = 0L)
  attempts$selected_primary <- "precision_retry"
  comp <- final$comp; colnames(comp$probabilities) <- levels(prep$data$Y)
  repaired <- .finish_fit(prep, candidate, comp, f$model_id, attempts,
    character(), character(), "ordinal", f$alternate)
  if (!isTRUE(repaired$valid) || !isTRUE(repaired$inference_valid))
    return(rejected("finish_fit_rejected", TRUE, polishes, checks))
  record$fit <- repaired; record$valid <- TRUE; record$reason <- ""
  record$values <- c(beta_B = unname(candidate["B_c"]), beta_A = unname(candidate["A_c"]))
  if (record$outcome == 6L) record$values <- c(record$values,
    d = unname(record$values["beta_A"] - record$values["beta_B"]))
  list(eligible = TRUE, accepted = TRUE, reason = "accepted_under_original_numerical_rules",
    initial_record = initial_record, record = record, polishes = polishes,
    checks = checks, protocol = mipr_protocol())
}

mipr_agreement_without_score <- function(theta, other, prep) {
  # Source primary is allowed to fail only its score gate. No other gate changes.
  first <- mipr_candidate(theta, prep, score_limit = Inf)
  second <- mipr_candidate(other, prep)
  if (!first$valid || !second$valid)
    return(list(valid = FALSE, reasons = c(first$reasons, second$reasons)))
  staff <- grep("^(B_c|A_c)", names(theta), value = TRUE)
  staff_difference <- if (length(staff)) max(abs(theta[staff] - other[staff])) else 0
  full_difference <- max(abs(theta - other)); ll_difference <- abs(first$comp$logLik - second$comp$logLik)
  se_ratio <- sqrt(diag(chol2inv(chol(second$comp$H))) / diag(chol2inv(chol(first$comp$H))))
  reasons <- c(if (staff_difference > 1e-4) "staffing_alternate_disagreement",
    if (full_difference > 0.01) "coefficient_growth_or_alternate_instability",
    if (ll_difference > 1e-6 * (1 + abs(first$comp$logLik))) "likelihood_alternate_disagreement",
    if (any(!is.finite(se_ratio)) || max(se_ratio) > 1.1) "standard_error_growth_under_tighter_tolerance")
  list(valid = !length(reasons), reasons = reasons, staff_difference = staff_difference,
    parameter_difference = full_difference, likelihood_difference = ll_difference,
    information_se_ratio = se_ratio)
}

mipr_overlay_registry <- function(registry, retries) {
  invalid <- which(!vapply(registry, function(x) isTRUE(x$valid), logical(1L)))
  indices <- vapply(retries, function(x) as.integer(x$index), integer(1L))
  if (anyNA(indices) || anyDuplicated(indices) || !setequal(indices, invalid))
    stop("Retry registry must cover every invalid source exactly once; subsets are forbidden.")
  output <- registry
  for (retry in retries) {
    index <- retry$index; source <- registry[[index]]
    if (!identical(source, retry$initial_record)) stop("Retry source record changed.")
    if (isTRUE(retry$accepted)) {
      if (!mipr_is_score_failure(source) || !isTRUE(retry$eligible) ||
          !isTRUE(retry$record$valid) || !isTRUE(retry$record$fit$valid) ||
          !isTRUE(retry$record$fit$inference_valid)) stop("Invalid accepted retry record.")
      for (nm in c("scenario", "seed", "chain", "iteration", "outcome", "n"))
        if (!identical(source[[nm]], retry$record[[nm]])) stop("Retry record identity changed.")
      output[[index]] <- retry$record
    } else if (!identical(source, retry$record)) stop("Rejected retry changed its source record.")
  }
  output
}
