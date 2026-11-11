# Specification 1.1, Section 8. Source this file; sourcing performs no analysis.
# Base/recommended R only. Numerical and chain failures remain in checkpoints.

mi_items <- paste0("Q", c(1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 15, 17,
                           21, 23, 24, 25, 26, 27))
mi_seeds <- c(20261111L, 20261112L)

mi_with_rng <- function(state = NULL, code) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  if (!is.null(state)) assign(".Random.seed", state, envir = .GlobalEnv)
  force(code)
}

mi_streams <- function(seed, m) {
  stopifnot(length(seed) == 1L, length(m) == 1L, m >= 1L, m == as.integer(m))
  mi_with_rng(code = {
    set.seed(seed)
    streams <- vector("list", m)
    streams[[1L]] <- get(".Random.seed", envir = .GlobalEnv)
    if (m > 1L) for (r in 2:m)
      streams[[r]] <- parallel::nextRNGStream(streams[[r - 1L]])
    streams
  })
}

mi_validate_data <- function(data) {
  if (!is.data.frame(data) || !all(c("Observation", "Site", mi_items) %in% names(data)))
    stop("MI input lacks required identifiers or questionnaire items.")
  if (anyNA(data$Observation) || anyDuplicated(data$Observation))
    stop("MI requires unique, nonmissing Observation identifiers.")
  site <- suppressWarnings(as.integer(as.character(data$Site)))
  if (anyNA(site) || !all(site %in% 1:8) || !setequal(site, 1:8))
    stop("MI requires known site and all eight sites.")
  for (item in mi_items) {
    v <- data[[item]]
    if (!is.numeric(v) || any(!is.na(v) & (!is.finite(v) | !v %in% 1:5)))
      stop(sprintf("MI requires decoded numeric values 1:5 or NA: %s.", item))
    if (anyNA(v) && sum(!is.na(v)) < 5L)
      stop(sprintf("Fewer than five originally observed donors for %s.", item))
  }
  invisible(TRUE)
}

mi_design <- function(completed, target) {
  others <- setdiff(mi_items, target)
  site <- as.integer(as.character(completed$Site))
  # One intercept, 18 numeric items centered at 3, and seven site indicators.
  # Centering changes neither the fitted values nor the specified posterior.
  X <- cbind(`(Intercept)` = 1, as.matrix(completed[others]) - 3,
             vapply(2:8, function(s) as.numeric(site == s), numeric(nrow(completed))))
  colnames(X) <- c("(Intercept)", others, paste0("Site", 2:8))
  storage.mode(X) <- "double"
  X
}

mi_regression_posterior <- function(X, y) {
  if (!is.matrix(X) || nrow(X) != length(y) || any(!is.finite(X)) ||
      any(!is.finite(y))) stop("Nonfinite or incompatible conditional regression input.")
  if (length(unique(y)) < 2L) stop("Conditional target is constant; recorded resolution required.")
  n <- nrow(X)
  q <- ncol(X)
  df <- n - q
  if (df <= 0L) stop("Conditional regression has nonpositive residual degrees of freedom.")
  decomposition <- qr(X, LAPACK = FALSE)
  if (decomposition$rank != q) stop("Conditional regression is rank deficient; amendment required.")
  estimate <- as.numeric(qr.coef(decomposition, y))
  residual <- as.numeric(qr.resid(decomposition, y))
  sse <- sum(residual^2)
  if (!is.finite(sse) || sse <= 0) stop("Conditional regression SSE is not positive.")
  R <- qr.R(decomposition)[seq_len(q), seq_len(q), drop = FALSE]
  # X[, pivot] = Q R. If z ~ N(0,I), R^{-1} z has covariance
  # R^{-1} R^{-T} = (R^T R)^{-1}; using R^{-T} z would be incorrect.
  sigma2 <- sse / rchisq(1L, df = df)
  if (!is.finite(sigma2) || sigma2 <= 0) stop("Invalid posterior variance draw.")
  displacement <- numeric(q)
  displacement[decomposition$pivot] <- sqrt(sigma2) * backsolve(R, rnorm(q))
  draw <- estimate + displacement
  if (any(!is.finite(draw))) stop("Nonfinite posterior coefficient draw.")
  list(estimate = estimate, draw = draw, sigma2 = sigma2, df = df,
       sse = sse, R = R, pivot = decomposition$pivot)
}

mi_choose_donors <- function(donor_prediction, recipient_prediction, donors = 5L) {
  if (donors != 5L || length(donor_prediction) < donors ||
      any(!is.finite(c(donor_prediction, recipient_prediction))))
    stop("Type-1 PMM requires five donors and finite predictions.")
  vapply(recipient_prediction, function(prediction) {
    distance <- abs(donor_prediction - prediction)
    # Random secondary ordering randomizes all exact boundary ties.
    closest <- order(distance, runif(length(distance)))[seq_len(donors)]
    closest[sample.int(donors, 1L)]
  }, integer(1L))
}

mi_update_item <- function(completed, original, target) {
  observed <- which(!is.na(original[[target]]))
  missing <- which(is.na(original[[target]]))
  if (!length(missing)) return(list(data = completed, donor_rows = integer()))
  if (length(observed) < 5L) stop("Insufficient originally observed donors.")
  X <- mi_design(completed, target)
  posterior <- mi_regression_posterior(X[observed, , drop = FALSE], original[[target]][observed])
  # Type 1: observed donor predictions use point estimates; recipient
  # predictions use the posterior draw. Imputed values never become donors.
  donor_prediction <- drop(X[observed, , drop = FALSE] %*% posterior$estimate)
  recipient_prediction <- drop(X[missing, , drop = FALSE] %*% posterior$draw)
  selected <- observed[mi_choose_donors(donor_prediction, recipient_prediction)]
  completed[[target]][missing] <- original[[target]][selected]
  list(data = completed, donor_rows = selected, recipient_rows = missing,
       df = posterior$df, sse = posterior$sse)
}

mi_initialize_chain <- function(original, stream, chain, seed) {
  mi_with_rng(stream, {
    completed <- original
    incomplete <- mi_items[vapply(original[mi_items], anyNA, logical(1L))]
    for (item in incomplete) {
      missing <- which(is.na(original[[item]]))
      observed <- original[[item]][!is.na(original[[item]])]
      completed[[item]][missing] <- observed[sample.int(length(observed), length(missing), replace = TRUE)]
    }
    list(chain = chain, seed = seed, initial_stream = stream,
         state = get(".Random.seed", envir = .GlobalEnv), data = completed,
         iteration = 0L, status = "ready", failure = NULL, trace = list(),
         donor_counts = setNames(lapply(incomplete, function(x) integer(nrow(original))), incomplete),
         terminal_donors = list())
  })
}

mi_advance_chain <- function(chain, original, iterations, checkpoint = NULL) {
  if (identical(chain$status, "failed")) return(chain)
  if (iterations < chain$iteration) stop("Cannot reduce completed chain iteration count.")
  incomplete <- mi_items[vapply(original[mi_items], anyNA, logical(1L))]
  if (chain$iteration == iterations) return(chain)
  mi_with_rng(chain$state, {
    for (iteration in seq.int(chain$iteration + 1L, iterations)) {
      target <- NA_character_
      failure <- tryCatch({
        for (target in incomplete) {
          update <- mi_update_item(chain$data, original, target)
          chain$data <- update$data
          chain$donor_counts[[target]] <- chain$donor_counts[[target]] +
            tabulate(update$donor_rows, nbins = nrow(original))
          chain$terminal_donors[[target]] <- data.frame(recipient_row = update$recipient_rows,
                                                       donor_row = update$donor_rows)
        }
        NULL
      }, error = function(e) conditionMessage(e))
      chain$state <- get(".Random.seed", envir = .GlobalEnv)
      if (!is.null(failure)) {
        chain$status <- "failed"
        chain$failure <- list(iteration = iteration, item = target, reason = failure)
        if (!is.null(checkpoint)) checkpoint(chain)
        return(chain)
      }
      chain$iteration <- iteration
      chain$status <- "complete"
      if (length(incomplete)) chain$trace[[iteration]] <- do.call(rbind, lapply(incomplete, function(item) {
        values <- chain$data[[item]][is.na(original[[item]])]
        data.frame(seed = chain$seed, chain = chain$chain, iteration = iteration,
                   item = item, n_imputed = length(values), mean = mean(values),
                   variance = if (length(values) > 1L) var(values) else NA_real_)
      }))
      if (!is.null(checkpoint) && (iteration %% 10L == 0L || iteration == iterations)) checkpoint(chain)
    }
    chain
  })
}

mi_new_checkpoint <- function(data) {
  mi_validate_data(data)
  list(schema = "pmm_specification_1.1_v1", original = data, seeds = mi_seeds,
       target_m = 50L, target_iterations = 50L, runs = list(), reviews = list(),
       fit_history = list(), fits = NULL, status = "generation_pending",
       created = format(Sys.time(), tz = "UTC", usetz = TRUE))
}

mi_save_checkpoint <- function(state, path) {
  # Keep the preceding checkpoint as a recovery copy. Neither copy is deleted.
  if (file.exists(path)) {
    if (!file.copy(path, paste0(path, ".previous"), overwrite = TRUE))
      stop("Could not preserve previous MI checkpoint.")
  }
  saveRDS(state, path, version = 3L)
  invisible(state)
}

mi_collect_trace <- function(state) {
  rows <- unlist(lapply(state$runs, function(run)
    lapply(run, function(chain) chain$trace)), recursive = FALSE)
  rows <- unlist(rows, recursive = FALSE)
  if (!length(rows)) return(data.frame(seed = integer(), chain = integer(),
    iteration = integer(), item = character(), n_imputed = integer(), mean = numeric(), variance = numeric()))
  do.call(rbind, rows)
}

mi_diagnostics <- function(state) {
  original <- state$original
  site <- as.integer(as.character(original$Site))
  distributions <- list()
  reuse <- list()
  correlations <- list()
  registry <- list()
  for (run in state$runs) for (chain in run) {
    registry[[length(registry) + 1L]] <- data.frame(seed = chain$seed, chain = chain$chain,
      iteration = chain$iteration, status = chain$status,
      failed_item = if (is.null(chain$failure)) NA_character_ else chain$failure$item,
      failed_iteration = if (is.null(chain$failure)) NA_integer_ else chain$failure$iteration,
      reason = if (is.null(chain$failure)) "" else chain$failure$reason)
    for (item in mi_items) for (s in 1:8) for (kind in c("observed", "imputed")) {
      selected <- site == s & (is.na(original[[item]]) == (kind == "imputed"))
      counts <- tabulate(as.integer(chain$data[[item]][selected]), nbins = 5L)
      distributions[[length(distributions) + 1L]] <- data.frame(seed = chain$seed,
        chain = chain$chain, iteration = chain$iteration, item = item, site = s,
        kind = kind, value = 1:5, n = counts, denominator = sum(selected))
    }
    for (item in names(chain$donor_counts)) {
      donors <- which(!is.na(original[[item]]))
      terminal <- chain$terminal_donors[[item]]
      final_counts <- if (is.null(terminal)) integer(nrow(original)) else
        tabulate(terminal$donor_row, nbins = nrow(original))
      reuse[[length(reuse) + 1L]] <- data.frame(seed = chain$seed, chain = chain$chain,
        iteration = chain$iteration, item = item, donor_row = donors,
        Observation = original$Observation[donors], donor_site = site[donors],
        uses_all_iterations = chain$donor_counts[[item]][donors],
        uses_terminal_iteration = final_counts[donors])
    }
    for (a in seq_along(mi_items)) if (a < length(mi_items))
      for (b in seq.int(a + 1L, length(mi_items))) {
        item_a <- mi_items[a]
        item_b <- mi_items[b]
        # Both completed-data correlations and correlations on rows involving
        # at least one imputed member are retained, with explicit denominators.
        for (scope in c("completed_all_rows", "at_least_one_imputed")) {
          selected <- if (scope == "completed_all_rows") rep(TRUE, nrow(original)) else
            is.na(original[[item_a]]) | is.na(original[[item_b]])
          va <- chain$data[[item_a]][selected]
          vb <- chain$data[[item_b]][selected]
          correlation <- if (length(va) >= 2L && sd(va) > 0 && sd(vb) > 0) cor(va, vb) else NA_real_
          correlations[[length(correlations) + 1L]] <- data.frame(seed = chain$seed,
            chain = chain$chain, iteration = chain$iteration, item_a = item_a,
            item_b = item_b, scope = scope, n = length(va), correlation = correlation)
        }
      }
  }
  bind <- function(rows) if (length(rows)) do.call(rbind, rows) else data.frame()
  list(trace = mi_collect_trace(state), distributions = bind(distributions),
       donor_reuse = bind(reuse), correlations = bind(correlations), chains = bind(registry))
}

mi_write_diagnostics <- function(state, outdir) {
  diagnostics <- mi_diagnostics(state)
  for (name in names(diagnostics))
    write.csv(diagnostics[[name]], file.path(outdir, paste0("mi_", name, ".csv")), row.names = FALSE)
  writeLines(c("PMM imputation diagnostics.",
    paste("Status:", state$status), paste("Chains per root seed:", state$target_m),
    paste("Target iterations:", state$target_iterations),
    paste("Monte Carlo status:", if (is.null(state$monte_carlo_status)) "not yet assessed" else state$monte_carlo_status),
    "Trace means/variances describe imputed entries after each complete iteration.",
    "Distributions, correlations, and terminal donor reuse describe saved chain states.",
    "Inspect drift and mixing before approving fits or extending the chain count.",
    "Between-imputation dispersion and Monte Carlo error are not sampling standard errors.",
    "No pooled tests, nominal intervals, or reduced-model fit gains are produced."),
    file.path(outdir, "mi_status.txt"))
  invisible(diagnostics)
}

mi_analysis_data <- function(original, completed, outcome) {
  if (!outcome %in% c(4L, 6L)) stop("MI supports Q6 and Q4 only.")
  selected <- !is.na(original[[paste0("Q", outcome)]])
  completed[selected, , drop = FALSE]
}

mi_default_fit <- function(data_completed, outcome, chain, seed) {
  if (!exists("derive_scores", mode = "function") || !exists("fit_ordinal", mode = "function"))
    stop("Source data_prep.R and model_engine.R before MI outcome fits.")
  prepared <- derive_scores(data_completed, outcome = outcome, management = "stored")
  fit_ordinal(prepared, outcome = "Y", predictors = c("L_c", "E_c", "R_c", "M_c", "B_c", "A_c"),
              site = "Site", model_id = sprintf("mi_Q%d_seed%d_chain%d", outcome, seed, chain))
}

mi_fit_checkpoint <- function(state, fit_function = mi_default_fit, checkpoint = NULL) {
  if (!identical(state$status, "diagnostics_accepted"))
    stop("A recorded acceptance of current chain diagnostics is required before fitting.")
  if (any(vapply(state$runs, function(run) length(run) != state$target_m, logical(1L))))
    stop("Required number of completed chains is unavailable.")
  results <- if (is.null(state$fits)) list() else state$fits
  for (run in state$runs) for (chain in run) for (outcome in c(6L, 4L)) {
    if (chain$status != "complete" || chain$iteration != state$target_iterations)
      stop("A required completed chain is unavailable; do not replace chains.")
    retained <- vapply(results, function(fit) fit$seed == chain$seed &&
      fit$chain == chain$chain && fit$outcome == outcome && fit$iteration == chain$iteration, logical(1L))
    if (any(retained)) next
    completed <- mi_analysis_data(state$original, chain$data, outcome)
    fit <- tryCatch(fit_function(completed, outcome, chain$chain, chain$seed),
                    error = function(e) list(valid = FALSE, status = "error", reason = conditionMessage(e)))
    coefficients <- fit$coefficients
    if (is.null(coefficients)) coefficients <- fit$coef
    valid <- isTRUE(fit$valid) && all(c("B_c", "A_c") %in% names(coefficients)) &&
      all(is.finite(coefficients[c("B_c", "A_c")]))
    values <- if (valid) c(beta_B = unname(coefficients["B_c"]), beta_A = unname(coefficients["A_c"])) else
      c(beta_B = NA_real_, beta_A = NA_real_)
    if (outcome == 6L) values <- c(values, d = unname(values["beta_A"] - values["beta_B"]))
    reason <- if (valid) "" else if (!is.null(fit$reason)) as.character(fit$reason) else
      if (!is.null(fit$status)) as.character(fit$status) else "Invalid full outcome fit."
    results[[length(results) + 1L]] <- list(seed = chain$seed, chain = chain$chain,
      iteration = chain$iteration, outcome = outcome, n = nrow(completed), valid = valid,
      reason = paste(reason, collapse = "; "), values = values, fit = fit)
    state$fits <- results
    if (!is.null(checkpoint)) checkpoint(state)
  }
  state$fits <- results
  state$status <- "fitted"
  state
}

mi_summarize_fits <- function(state) {
  if (is.null(state$fits)) stop("No imputed-data outcome fits are available.")
  estimates <- list()
  covariance <- list()
  summary <- list()
  for (seed in state$seeds) for (outcome in c(6L, 4L)) {
    selected <- Filter(function(fit) fit$seed == seed && fit$outcome == outcome, state$fits)
    n_valid <- sum(vapply(selected, function(fit) fit$valid, logical(1L)))
    valid <- length(selected) == state$target_m && n_valid == state$target_m
    summary[[length(summary) + 1L]] <- data.frame(seed = seed, outcome = outcome,
      required = state$target_m, attempted = length(selected), valid = n_valid,
      failed = length(selected) - n_valid, status = if (valid) "available" else "unavailable")
    if (!valid) next
    values <- do.call(rbind, lapply(selected, function(fit) fit$values))
    average <- colMeans(values)
    between <- cov(values)
    estimates[[length(estimates) + 1L]] <- data.frame(seed = seed, outcome = outcome,
      component = names(average), m = nrow(values), estimate = unname(average),
      exponentiated = exp(unname(average)),
      exponentiated_quantity = ifelse(names(average) == "d", "ratio_of_odds_ratios", "odds_ratio"),
      between_imputation_variance = unname(diag(between)),
      monte_carlo_error = sqrt(unname(diag(between)) / nrow(values)))
    indices <- expand.grid(component_row = colnames(between), component_column = colnames(between),
                           stringsAsFactors = FALSE)
    covariance[[length(covariance) + 1L]] <- data.frame(seed = seed, outcome = outcome,
      indices, between_imputation_covariance = as.vector(between))
  }
  bind <- function(rows) if (length(rows)) do.call(rbind, rows) else data.frame()
  estimates <- bind(estimates)
  # Failure in either seed makes the overall sensitivity unavailable for that
  # outcome. Do not report a retained successful-seed aggregate as its result.
  summary <- bind(summary)
  available_outcomes <- c(6L, 4L)[vapply(c(6L, 4L), function(outcome)
    all(summary$status[summary$outcome == outcome] == "available"), logical(1L))]
  if (nrow(estimates)) estimates <- estimates[estimates$outcome %in% available_outcomes, , drop = FALSE]
  covariance <- bind(covariance)
  if (nrow(covariance)) covariance <- covariance[covariance$outcome %in% available_outcomes, , drop = FALSE]
  checks <- list()
  if (nrow(estimates)) for (outcome in available_outcomes) {
    first <- estimates[estimates$seed == state$seeds[1L] & estimates$outcome == outcome, , drop = FALSE]
    second <- estimates[estimates$seed == state$seeds[2L] & estimates$outcome == outcome, , drop = FALSE]
    second <- second[match(first$component, second$component), , drop = FALSE]
    delta <- abs(first$estimate - second$estimate)
    max_mc <- pmax(first$monte_carlo_error, second$monte_carlo_error)
    checks[[length(checks) + 1L]] <- data.frame(outcome = outcome, component = first$component,
      absolute_seed_difference = delta, maximum_monte_carlo_error = max_mc,
      exceeds_tolerance = delta > 0.02 | max_mc > 0.01)
  }
  checks <- bind(checks)
  if (!nrow(estimates)) estimates <- data.frame(seed = integer(), outcome = integer(),
    component = character(), m = integer(), estimate = numeric(), exponentiated = numeric(),
    exponentiated_quantity = character(), between_imputation_variance = numeric(), monte_carlo_error = numeric())
  if (!nrow(covariance)) covariance <- data.frame(seed = integer(), outcome = integer(),
    component_row = character(), component_column = character(), between_imputation_covariance = numeric())
  if (!nrow(checks)) checks <- data.frame(outcome = integer(), component = character(),
    absolute_seed_difference = numeric(), maximum_monte_carlo_error = numeric(), exceeds_tolerance = logical())
  list(estimates = estimates, covariance = covariance, availability = summary,
       monte_carlo_checks = checks,
       extension_required = nrow(checks) > 0L && any(checks$exceeds_tolerance),
       at_chain_cap = state$target_m == 200L)
}

mi_write_fits <- function(state, outdir) {
  summaries <- mi_summarize_fits(state)
  for (name in c("estimates", "covariance", "availability", "monte_carlo_checks"))
    write.csv(summaries[[name]], file.path(outdir, paste0("mi_", name, ".csv")), row.names = FALSE)
  registry <- do.call(rbind, lapply(state$fits, function(fit) data.frame(seed = fit$seed,
    chain = fit$chain, iteration = fit$iteration, outcome = fit$outcome, n = fit$n,
    valid = fit$valid, reason = fit$reason)))
  write.csv(registry, file.path(outdir, "mi_model_registry.csv"), row.names = FALSE)
  # Complete solver attempts and warnings are retained in state$fits[[i]]$fit.
  saveRDS(summaries, file.path(outdir, "mi_summary.rds"), version = 3L)
  invisible(summaries)
}

mi_review <- function(review, allowed) {
  if (!is.list(review) || !all(c("reviewer", "notes", "decision") %in% names(review)) ||
      any(!vapply(review[c("reviewer", "notes", "decision")], function(value)
        is.character(value) && length(value) == 1L && !is.na(value) && nzchar(trimws(value)), logical(1L))) ||
      !review$decision %in% allowed)
    stop(paste("Supply reviewer, notes, and a decision from:", paste(allowed, collapse = ", ")))
  review$recorded <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  review
}

# Initial call: run_mi(data, outdir). Inspect saved diagnostics, then call again
# with diagnostic_review=list(reviewer="...", notes="...", decision="accept").
# Iteration extensions require decision="extend_iterations" and iterations=100
# then 200. Chain extensions require decision="extend_chains" and m=100 then
# 200, after a fitted Monte Carlo check exceeds the stated tolerance. They stop
# for inspection of the new chains before fitting. The same output directory
# resumes the checkpoint; different inputs cannot reuse it.
run_mi <- function(data, outdir, m = 50L, iterations = 50L,
                   fit_function = NULL, diagnostic_review = NULL) {
  mi_validate_data(data)
  if (length(m) != 1L || !m %in% c(50L, 100L, 200L) || length(iterations) != 1L ||
      !iterations %in% c(50L, 100L, 200L)) stop("Planned chain and iteration counts are 50, 100, or 200.")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(outdir, "mi_checkpoint.rds")
  existing <- file.exists(path)
  state <- if (existing) readRDS(path) else mi_new_checkpoint(data)
  if (!identical(state$schema, "pmm_specification_1.1_v1") || !identical(state$original, data))
    stop("Checkpoint schema or original input differs; use a distinct output directory.")
  if (!existing && (m != 50L || iterations != 50L)) stop("Begin with 50 chains and 50 iterations.")
  if (identical(state$status, "unavailable")) stop("Sensitivity was marked unavailable; amendment required.")
  review <- NULL
  extension <- m != state$target_m || iterations != state$target_iterations
  if (extension) {
    if (m < state$target_m || iterations < state$target_iterations)
      stop("Cannot discard completed chains or iterations.")
    if (m != state$target_m && iterations != state$target_iterations)
      stop("Review iteration and chain-count extensions separately.")
    if (iterations != state$target_iterations) {
      review <- mi_review(diagnostic_review, "extend_iterations")
      if (iterations != 2L * state$target_iterations || state$target_iterations >= 200L)
        stop("Iteration extension must proceed from 50 to 100 to 200.")
    } else {
      review <- mi_review(diagnostic_review, "extend_chains")
      if (m != 2L * state$target_m || state$target_m >= 200L)
        stop("Chain extension must proceed from 50 to 100 to 200.")
      summaries <- mi_summarize_fits(state)
      if (!summaries$extension_required) stop("No specified Monte Carlo criterion calls for more chains.")
    }
    if (!is.null(state$fits)) state$fit_history[[length(state$fit_history) + 1L]] <- state$fits
    # Additional chains leave the earlier completed data and full-fit outcomes
    # intact, including failures. More iterations change every completed data
    # set, so only that extension requires new fits on all terminal data sets.
    if (iterations != state$target_iterations) state$fits <- NULL
    state$monte_carlo_status <- NULL
    state$target_m <- as.integer(m)
    state$target_iterations <- as.integer(iterations)
    state$status <- "generation_pending"
    state$reviews[[length(state$reviews) + 1L]] <- c(review, list(m = m, iterations = iterations))
  }
  if (state$status %in% c("generation_pending", "generating")) {
    state$status <- "generating"
    for (seed in state$seeds) {
      key <- as.character(seed)
      streams <- mi_streams(seed, state$target_m)
      if (is.null(state$runs[[key]])) state$runs[[key]] <- list()
      for (r in seq_len(state$target_m)) {
        if (length(state$runs[[key]]) < r)
          state$runs[[key]][[r]] <- mi_initialize_chain(data, streams[[r]], r, seed)
        # Lexical state updates preserve intermediate states without consuming RNG.
        save_chain <- function(value) {
          state$runs[[key]][[r]] <<- value
          mi_save_checkpoint(state, path)
        }
        state$runs[[key]][[r]] <- mi_advance_chain(state$runs[[key]][[r]], data,
                                                   state$target_iterations, save_chain)
        mi_save_checkpoint(state, path)
        if (state$runs[[key]][[r]]$status == "failed") {
          state$status <- "conditional_model_failure"
          mi_save_checkpoint(state, path)
          mi_write_diagnostics(state, outdir)
          return(invisible(state))
        }
      }
    }
    state$status <- "diagnostics_pending"
    mi_save_checkpoint(state, path)
    mi_write_diagnostics(state, outdir)
    # Every newly generated terminal collection requires its own inspection.
    return(invisible(state))
  }
  if (state$status == "conditional_model_failure")
    stop("A conditional model failed. Inspect the saved failure; recorded resolution/amendment required.")
  if (!is.null(diagnostic_review)) {
    review <- mi_review(diagnostic_review, c("accept", "unavailable"))
    if (state$status != "diagnostics_pending") stop("Current state is not awaiting diagnostic review.")
    if (review$decision == "unavailable" && state$target_iterations != 200L)
      stop("Unresolved drift must first extend all chains to 100 and then 200 iterations.")
    state$reviews[[length(state$reviews) + 1L]] <- c(review,
      list(m = state$target_m, iterations = state$target_iterations))
    state$status <- if (review$decision == "accept") "diagnostics_accepted" else "unavailable"
    mi_save_checkpoint(state, path)
  }
  if (state$status == "diagnostics_accepted") {
    if (is.null(fit_function)) fit_function <- mi_default_fit
    state <- mi_fit_checkpoint(state, fit_function, function(value) mi_save_checkpoint(value, path))
    mi_save_checkpoint(state, path)
  }
  if (state$status == "fitted") {
    summaries <- mi_write_fits(state, outdir)
    state$monte_carlo_status <- if (!nrow(summaries$estimates)) "unavailable_no_valid_outcome" else
      if (!summaries$extension_required) "within_tolerances_for_available_outcomes" else
      if (summaries$at_chain_cap) "instability_at_200_chain_cap" else "extension_required"
    mi_save_checkpoint(state, path)
  }
  mi_write_diagnostics(state, outdir)
  invisible(state)
}
