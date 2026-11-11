# Synthetic verification only. This file defines tests without executing them.
# After the implementation review: source imputation.R, source this file, then
# run_imputation_tests(). No study input or study outcome model is used here.

run_imputation_tests <- function() {
  require_function <- function(name) {
    if (!exists(name, mode = "function")) stop("Source imputation.R first: ", name)
  }
  require_function("mi_update_item")
  expect_error <- function(expression, pattern) {
    error <- tryCatch({ force(expression); NULL }, error = function(e) conditionMessage(e))
    if (is.null(error) || !grepl(pattern, error)) stop("Expected error matching: ", pattern)
    invisible(TRUE)
  }
  checks <- character()
  passed <- function(name) checks <<- c(checks, name)
  synthetic <- mi_with_rng(mi_streams(901L, 1L)[[1L]], {
    n <- 240L
    data <- data.frame(Observation = as.character(seq_len(n)), Site = factor(rep(1:8, each = n / 8)))
    for (item in mi_items) data[[item]] <- sample.int(5L, n, replace = TRUE)
    for (i in seq_along(mi_items)) data[[mi_items[i]]][c(i, i + 35L, i + 100L)] <- NA_integer_
    data
  })
  mi_validate_data(synthetic)
  streams <- mi_streams(20261111L, 5L)
  stopifnot(identical(streams, mi_streams(20261111L, 5L)),
            identical(streams, mi_streams(20261111L, 9L)[seq_len(5L)]),
            !identical(streams[[1L]], streams[[2L]]),
            !identical(streams[[1L]], mi_streams(20261112L, 1L)[[1L]]))
  passed("deterministic distinct streams and chain-extension prefix")

  # Helpers restore the caller's random generator and state.
  original_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) original_seed <- get(".Random.seed", envir = .GlobalEnv)
  invisible(mi_streams(20261111L, 2L))
  stopifnot(identical(RNGkind(), original_kind),
            identical(had_seed, exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)))
  if (had_seed) stopifnot(identical(original_seed, get(".Random.seed", envir = .GlobalEnv)))
  passed("caller RNG state is preserved")

  initial <- mi_initialize_chain(synthetic, streams[[1L]], 1L, 20261111L)
  direct <- mi_advance_chain(initial, synthetic, 4L)
  resumed <- mi_advance_chain(mi_advance_chain(initial, synthetic, 2L), synthetic, 4L)
  repeated <- mi_advance_chain(mi_initialize_chain(synthetic, streams[[1L]], 1L, 20261111L), synthetic, 4L)
  stopifnot(direct$status == "complete", identical(direct, resumed), identical(direct, repeated))
  for (item in mi_items) {
    observed <- !is.na(synthetic[[item]])
    stopifnot(identical(direct$data[[item]][observed], synthetic[[item]][observed]),
              all(direct$data[[item]] %in% 1:5),
              all(direct$data[[item]][!observed] %in% synthetic[[item]][observed]))
    mapping <- direct$terminal_donors[[item]]
    stopifnot(all(observed[mapping$donor_row]), all(!observed[mapping$recipient_row]),
              all(direct$data[[item]][mapping$recipient_row] == synthetic[[item]][mapping$donor_row]),
              sum(direct$donor_counts[[item]]) == 4L * sum(!observed))
  }
  stopifnot(identical(direct$data$Site, synthetic$Site),
            identical(direct$data$Observation, synthetic$Observation))
  passed("observed values unchanged, observed donors only, support preserved, exact resume")

  # Explicitly reconstruct the Type-1 operation with the identical RNG state.
  target <- "Q6"
  observed <- which(!is.na(synthetic[[target]]))
  missing <- which(is.na(synthetic[[target]]))
  X <- mi_design(initial$data, target)
  stopifnot(ncol(X) == 26L, sum(colnames(X) == "(Intercept)") == 1L,
            identical(colnames(X)[20:26], paste0("Site", 2:8)))
  expected_donors <- mi_with_rng(streams[[2L]], {
    posterior <- mi_regression_posterior(X[observed, , drop = FALSE], synthetic[[target]][observed])
    observed[mi_choose_donors(drop(X[observed, , drop = FALSE] %*% posterior$estimate),
                             drop(X[missing, , drop = FALSE] %*% posterior$draw))]
  })
  actual <- mi_with_rng(streams[[2L]], mi_update_item(initial$data, synthetic, target))
  stopifnot(identical(actual$donor_rows, expected_donors))
  passed("Type-1 donor point prediction and recipient posterior prediction")

  # All six donors tie at the boundary; randomized candidate selection must
  # permit every donor, including the donor outside a fixed first-five slice.
  tied <- mi_with_rng(streams[[3L]], mi_choose_donors(rep(0, 6L), rep(0, 600L)))
  stopifnot(setequal(tied, 1:6))
  nearest <- mi_with_rng(streams[[3L]], mi_choose_donors(c(0:4, 100, 200), rep(0, 200L)))
  stopifnot(all(nearest %in% 1:5))
  passed("five nearest donors and random boundary ties")

  # Nonorthogonal design detects the inverse-triangular orientation error.
  test_x <- cbind(1, seq(-2, 3, length.out = 40L),
                   sin(seq_len(40L)) + seq(-1, 2, length.out = 40L))
  test_y <- 1 + test_x[, 2] - 0.6 * test_x[, 3] + cos(seq_len(40L) * 1.7)
  posterior <- mi_with_rng(streams[[4L]], mi_regression_posterior(test_x, test_y))
  inverse <- solve(crossprod(test_x))
  inverse_r <- backsolve(posterior$R, diag(ncol(test_x)))
  from_qr <- matrix(0, ncol(test_x), ncol(test_x))
  from_qr[posterior$pivot, posterior$pivot] <- tcrossprod(inverse_r)
  stopifnot(max(abs(inverse - from_qr)) < 1e-10,
            max(abs(posterior$estimate - solve(crossprod(test_x), crossprod(test_x, test_y)))) < 1e-10,
            posterior$df == nrow(test_x) - ncol(test_x),
            abs(posterior$sse - sum((test_y - test_x %*% posterior$estimate)^2)) < 1e-10)
  standardized <- mi_with_rng(streams[[4L]], replicate(4000L, {
    draw <- mi_regression_posterior(test_x, test_y)
    as.numeric(chol(crossprod(test_x)) %*% ((draw$draw - draw$estimate) / sqrt(draw$sigma2)))
  }))
  stopifnot(max(abs(rowMeans(standardized))) < 0.06,
            max(abs(cov(t(standardized)) - diag(3L))) < 0.09)
  passed("posterior mean, SSE, n-minus-q degrees of freedom, and covariance orientation")

  # Verify the variance draw directly from the same stream and chi-square law.
  sigma_expected <- mi_with_rng(streams[[5L]], posterior$sse / rchisq(1L, posterior$df))
  sigma_actual <- mi_with_rng(streams[[5L]], mi_regression_posterior(test_x, test_y)$sigma2)
  stopifnot(identical(sigma_expected, sigma_actual))
  expect_error(mi_regression_posterior(cbind(test_x, test_x[, 2L]), test_y), "rank deficient")
  expect_error(mi_regression_posterior(test_x, rep(3, nrow(test_x))), "constant")
  expect_error(mi_regression_posterior(diag(3), 1:3), "nonpositive")
  passed("specified variance posterior and mandatory conditional-model failures")

  for (outcome in c(6L, 4L)) {
    analysis <- mi_analysis_data(synthetic, direct$data, outcome)
    expected_ids <- synthetic$Observation[!is.na(synthetic[[paste0("Q", outcome)]])]
    stopifnot(identical(analysis$Observation, expected_ids),
              !anyNA(analysis[[paste0("Q", outcome)]]))
  }
  passed("originally missing outcomes excluded after auxiliary imputation")

  # Fit callbacks are synthetic constants. No ordinal outcome model is fit.
  state <- mi_new_checkpoint(synthetic)
  state$target_m <- 2L
  state$target_iterations <- 4L
  state$status <- "diagnostics_accepted"
  for (seed in mi_seeds) {
    a <- direct
    a$seed <- seed
    b <- direct
    b$seed <- seed
    b$chain <- 2L
    state$runs[[as.character(seed)]] <- list(a, b)
  }
  fake_fit <- function(data_completed, outcome, chain, seed) {
    expected_ids <- synthetic$Observation[!is.na(synthetic[[paste0("Q", outcome)]])]
    stopifnot(identical(data_completed$Observation, expected_ids))
    list(valid = TRUE, coefficients = c(B_c = chain * 0.02, A_c = chain * 0.03), status = "valid")
  }
  fitted <- mi_fit_checkpoint(state, fake_fit)
  summary <- mi_summarize_fits(fitted)
  stopifnot(nrow(summary$estimates) == 10L,
            all(summary$availability$status == "available"),
            abs(summary$estimates$estimate[summary$estimates$component == "d"][1] - 0.015) < 1e-12)
  fake_failure <- function(data_completed, outcome, chain, seed) {
    if (outcome == 6L && chain == 1L && seed == mi_seeds[1L])
      return(list(valid = FALSE, status = "synthetic_failure"))
    fake_fit(data_completed, outcome, chain, seed)
  }
  failed <- mi_fit_checkpoint(state, fake_failure)
  failure_summary <- mi_summarize_fits(failed)
  stopifnot(!any(failure_summary$estimates$outcome == 6L),
            all(failure_summary$estimates$outcome == 4L),
            sum(failure_summary$availability$failed) == 1L)
  passed("point averaging, Q6 difference, and no successful-fit subset pooling")

  # A checkpointed failed fit must remain failed; resumption cannot replace it.
  failed$status <- "diagnostics_accepted"
  retained <- mi_fit_checkpoint(failed, function(...) stop("Must retain completed fit attempts."))
  stopifnot(identical(failed$fits, retained$fits))
  expect_error(mi_fit_checkpoint(mi_new_checkpoint(synthetic), fake_fit), "recorded acceptance")
  passed("fit resumption retains failures and requires diagnostic review")

  diagnostics <- mi_diagnostics(state)
  stopifnot(nrow(diagnostics$trace) == 2L * 2L * 4L * length(mi_items),
            nrow(diagnostics$chains) == 4L,
            all(diagnostics$distributions$n >= 0),
            all(diagnostics$donor_reuse$uses_all_iterations >= 0))
  passed("diagnostic trace and terminal output dimensions")
  data.frame(check = checks, passed = TRUE, stringsAsFactors = FALSE)
}
