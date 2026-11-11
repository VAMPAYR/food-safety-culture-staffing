# Synthetic verification for specification 1.4. Defining these tests runs no R
# analysis. Execute only after the independent specification and code reviews:
# source("imputation_extensions.R"); source("test_imputation_extensions.R")
# run_imputation_extension_tests(script_dir = "C:/absolute/path/to/reproducibility/R")
# No study records or ordinal outcome fits are used by this test suite.

# Separate, small manifest test for the public-clone packaging fix. This does
# not fit models, generate imputations, or change the main suite's API.
run_imputation_manifest_tests <- function(script_dir, fixture_dir) {
  script_dir <- normalizePath(script_dir,winslash="/",mustWork=TRUE)
  if(dir.exists(fixture_dir) && length(list.files(fixture_dir,all.files=TRUE,no..=TRUE))) stop("Manifest fixture directory must be new or empty.")
  dir.create(fixture_dir,recursive=TRUE,showWarnings=FALSE)
  isolated <- new.env(parent=globalenv())
  source(file.path(script_dir,"run_imputation_extensions.R"),local=isolated)
  required <- c("data_prep.R","model_engine.R","imputation_extensions.R","run_imputation_extensions.R","test_imputation_extensions.R")
  public_r <- file.path(fixture_dir,"public","R"); dir.create(public_r,recursive=TRUE)
  for(name in required) writeLines("# Manifest-only fixture; never executed.",file.path(public_r,name))
  public_paths <- isolated$mie_source_paths(public_r)
  stopifnot(length(public_paths)==5L,setequal(basename(public_paths),required),all(file.exists(public_paths)))
  internal_spec <- file.path(fixture_dir,"public","specification","mathematical_specification","main.tex")
  dir.create(dirname(internal_spec),recursive=TRUE); writeLines("Internal manifest test fixture.",internal_spec)
  local_paths <- isolated$mie_source_paths(public_r)
  stopifnot(length(local_paths)==6L,identical(local_paths[seq_along(public_paths)],public_paths),identical(tail(local_paths,1L),internal_spec))
  incomplete_r <- file.path(fixture_dir,"incomplete","R"); dir.create(incomplete_r,recursive=TRUE)
  for(name in required[-1L]) writeLines("# Manifest-only fixture.",file.path(incomplete_r,name))
  message <- tryCatch({isolated$mie_source_paths(incomplete_r);NULL},error=function(e)conditionMessage(e))
  stopifnot(!is.null(message),grepl("required reviewed R source",message))
  data.frame(test=c("public source manifest excludes absent internal specification","local source manifest includes existing internal specification","missing executable source remains an error"),passed=TRUE,detail="Manifest fixture passed; no analysis ran.",stringsAsFactors=FALSE)
}

run_imputation_extension_tests <- function(script_dir) {
  if (missing(script_dir) || length(script_dir) != 1L || !dir.exists(script_dir))
    stop("Supply the absolute R source directory as script_dir.")
  script_dir <- normalizePath(script_dir, winslash = "/", mustWork = TRUE)
  if (!all(file.exists(file.path(script_dir, c("data_prep.R", "model_engine.R", "imputation_extensions.R")))))
    stop("The supplied directory lacks the required R source files.")
  needed <- c("mie_variables", "mie_encode", "mie_validate", "mie_design",
              "mie_posterior", "mie_weights", "mie_draw_donors", "mie_streams",
              "mie_with_rng", "mie_initialize", "mie_advance", "mie_analysis_data",
              "mie_scenarios", "mie_update", "mie_summarize_fits", "mie_tasks",
              "mie_parallel", "mie_generation_worker", "mie_write_diagnostics",
              "mie_validate_diagnostics", "mie_require_review", "mie_load_saved_fit", "mie_fit_tasks")
  for (name in needed)
    if (!exists(name, mode = "function")) stop("Source imputation_extensions.R first: ", name)
  checks <- character()
  passed <- function(name) checks <<- c(checks, name)
  expect_error <- function(expression, pattern = NULL) {
    message <- tryCatch({ force(expression); NULL }, error = function(e) conditionMessage(e))
    if (is.null(message) || (!is.null(pattern) && !grepl(pattern, message)))
      stop("Expected an error", if (!is.null(pattern)) paste("matching", pattern))
    invisible(TRUE)
  }
  close <- function(actual, expected, tolerance = 1e-10) {
    if (!identical(dim(actual), dim(expected)) || length(actual) != length(expected) ||
        any(!is.finite(c(actual, expected))) || max(abs(actual - expected)) > tolerance)
      stop("Independent numerical comparison failed.")
    invisible(TRUE)
  }

  expected_variables <- c(paste0("Q", c(1, 2, 3, 4, 6, 7, 8, 9, 11, 12,
                                       13, 15, 17, 21, 23, 24, 25, 26, 27)),
                          "D1", "D2", "D5")
  stopifnot(identical(mie_variables(), expected_variables))
  passed("fixed 22-variable membership and update order")

  streams <- mie_streams(20261111L, 6L)
  stopifnot(identical(streams, mie_streams(20261111L, 6L)),
            identical(streams, mie_streams(20261111L, 10L)[seq_len(6L)]),
            !identical(streams[[1L]], streams[[2L]]),
            !identical(streams[[1L]], mie_streams(20261112L, 1L)[[1L]]))
  original_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) original_seed <- get(".Random.seed", envir = .GlobalEnv)
  invisible(mie_with_rng(streams[[1L]], rnorm(7L)))
  stopifnot(identical(original_kind, RNGkind()),
            identical(had_seed, exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)))
  if (had_seed) stopifnot(identical(original_seed, get(".Random.seed", envir = .GlobalEnv)))
  passed("distinct deterministic streams, extension prefixes, and caller RNG preservation")

  # A nonorthogonal design exposes the inverse-Cholesky orientation error.
  # This calculation uses a direct inverse only as an independent small test.
  X <- cbind(1, seq(-0.5, 1.5, length.out = 24L),
             cos(seq_len(24L) * 0.7) + seq(0, 0.6, length.out = 24L))
  z <- 0.4 + 0.08 * X[, 2L] - 0.05 * X[, 3L] + 0.03 * sin(seq_len(24L))
  p <- ncol(X)
  m0 <- c(0.5, rep(0, p - 1L))
  precision <- crossprod(X) + diag(0.01, p)
  inverse <- solve(precision)
  expected_mean <- as.numeric(inverse %*% (crossprod(X, z) + 0.01 * m0))
  expected_shape <- 2 + nrow(X) / 2
  expected_scale <- 0.04 + 0.5 *
    (sum((z - as.numeric(X %*% expected_mean))^2) +
     0.01 * sum((expected_mean - m0)^2))
  actual <- mie_with_rng(streams[[1L]], mie_posterior(X, z))
  close(actual$precision, precision)
  close(as.numeric(actual$mean), expected_mean)
  close(actual$shape, expected_shape)
  close(actual$scale, expected_scale)
  inverse_chol <- backsolve(actual$chol, diag(p))
  close(tcrossprod(inverse_chol), inverse)
  stopifnot(max(abs(crossprod(inverse_chol) - inverse)) > 1e-4)

  expected_draw <- mie_with_rng(streams[[1L]], {
    v <- expected_scale / rgamma(1L, shape = expected_shape, rate = 1)
    list(sigma2 = v, draw = expected_mean + sqrt(v) *
           as.numeric(backsolve(chol(precision), rnorm(p))))
  })
  close(actual$sigma2, expected_draw$sigma2)
  close(as.numeric(actual$draw), expected_draw$draw)
  passed("independent proper posterior mean, precision, shape, scale, and draw convention")

  standardized <- mie_with_rng(streams[[2L]], replicate(4000L, {
    posterior <- mie_posterior(X, z)
    as.numeric(chol(precision) %*%
                 ((posterior$draw - expected_mean) / sqrt(posterior$sigma2)))
  }))
  stopifnot(max(abs(rowMeans(standardized))) < 0.07,
            max(abs(cov(t(standardized)) - diag(p))) < 0.10)
  passed("posterior coefficient draw has the specified covariance orientation")

  # Exact duplicate and all-zero columns must remain in the proper posterior.
  singular_X <- cbind(X, duplicate = X[, 2L], absent_category = 0)
  singular <- mie_with_rng(streams[[3L]], mie_posterior(singular_X, z))
  stopifnot(length(singular$mean) == ncol(singular_X),
            length(singular$draw) == ncol(singular_X),
            all(is.finite(c(singular$mean, singular$draw, singular$sigma2))),
            min(eigen(singular$precision, symmetric = TRUE, only.values = TRUE)$values) > 0)
  close(singular$precision - crossprod(singular_X), diag(0.01, ncol(singular_X)))
  close(singular$mean[ncol(singular_X)], 0)
  close(singular$shape, 2 + nrow(singular_X) / 2)
  # Binary and five-level target regressions share the same constants after
  # the prescribed target scaling. No data-dependent prior is permitted.
  binary_z <- rep(c(0, 1), length.out = nrow(X))
  binary <- mie_with_rng(streams[[3L]], mie_posterior(X, binary_z))
  close(binary$precision - crossprod(X), diag(0.01, p))
  close(binary$shape, expected_shape)
  close(binary$scale - 0.5 * (sum((binary_z - X %*% binary$mean)^2) +
                              0.01 * sum((binary$mean - m0)^2)), 0.04)
  bad_X <- X
  bad_X[1L, 2L] <- Inf
  expect_error(mie_posterior(bad_X, z))
  expect_error(mie_posterior(X, z[-1L]))
  passed("proper prior retains collinear and zero columns with unchanged constants")

  values <- 1:5
  uniform <- mie_weights(values, 0)
  up <- mie_weights(values, log(2))
  down <- mie_weights(values, -log(2))
  close(as.numeric(uniform), rep(0.2, 5L))
  close(as.numeric(up[-1L] / up[-5L]), rep(2, 4L))
  close(as.numeric(down[-1L] / down[-5L]), rep(0.5, 4L))
  stopifnot(sum(values * up) > sum(values * uniform),
            sum(values * down) < sum(values * uniform),
            sum(up[values >= 4L]) > sum(uniform[values >= 4L]),
            sum(down[values >= 4L]) < sum(uniform[values >= 4L]))
  close(as.numeric(mie_weights(rep(4L, 5L), log(2))), rep(0.2, 5L))
  close(as.numeric(mie_weights(rep(4L, 5L), -log(2))), rep(0.2, 5L))
  extreme <- mie_weights(values, 1000)
  stopifnot(all(is.finite(extreme)), all(extreme >= 0), abs(sum(extreme) - 1) < 1e-12,
            extreme[5L] > 0.999999)
  duplicated <- c(1, 2, 2, 4, 5)
  expected_weights <- 2^(duplicated - 3)
  expected_weights <- expected_weights / sum(expected_weights)
  close(as.numeric(mie_weights(duplicated, log(2))), expected_weights)
  passed("uniform baseline, adjacent-category tilt ratios, local support, and stable weights")

  nearest <- mie_with_rng(streams[[3L]], mie_draw_donors(
    c(0:4, 100, 200), rep(0, 200L), c(1:5, 1, 2), delta = 0))
  stopifnot(all(nearest$selected %in% 1:5),
            identical(dim(nearest$sets), c(200L, 5L)),
            all(apply(nearest$sets, 1L, function(x) setequal(x, 1:5))))
  tied <- mie_with_rng(streams[[3L]], mie_draw_donors(
    rep(0, 6L), rep(0, 600L), rep(3L, 6L), delta = log(2)))
  stopifnot(setequal(tied$selected, 1:6),
            all(tied$details$unique_categories == 1L),
            all(tied$details$single_category),
            all(abs(tied$details$conditional_mean_shift) < 1e-12),
            all(abs(tied$details$agreement_pp_shift) < 1e-12))
  tilted <- mie_with_rng(streams[[3L]], mie_draw_donors(
    0:4, rep(0, 5L), 1:5, delta = log(2)))
  close(tilted$details$conditional_mean_shift,
        rep(sum((1:5) * up) - 3, 5L))
  close(tilted$details$agreement_pp_shift,
        rep(100 * (sum(up[4:5]) - 0.4), 5L))
  stopifnot(all(tilted$details$unique_categories == 5L),
            all(as.matrix(tilted$details[, paste0("n", 1:5), drop = FALSE]) == 1))
  expect_error(mie_draw_donors(1:4, 1, 1:4))
  passed("five nearest observed donors, random boundary ties, and interpretable tilt diagnostics")

  # The fixture has all eight sites, genuine auxiliary gaps, and an entirely
  # absent predictor category. It deliberately contains exact collinearities.
  raw <- mie_with_rng(streams[[4L]], {
    n <- 200L
    d <- data.frame(Observation = as.character(seq_len(n)),
                    Site = factor(rep(1:8, each = 25L), levels = 1:8),
                    Language = factor(rep(c("English", "Spanish"), length.out = n),
                                      levels = c("English", "Spanish")))
    for (item in expected_variables[grepl("^Q", expected_variables)])
      d[[item]] <- sample.int(5L, n, replace = TRUE)
    tenure <- c("Less than 1 year", "1-5 years", "6-10 years", "11-15 years",
                "Greater than 15 years")
    d$D1 <- factor(tenure[sample.int(5L, n, replace = TRUE)], levels = tenure)
    d$D5 <- factor(tenure[sample.int(5L, n, replace = TRUE)], levels = tenure)
    d$D2 <- factor(rep(c("false", "true"), length.out = n), levels = c("false", "true"))
    d$Q27 <- rep(1L, n)
    d$D8 <- rep("never_an_imputation_predictor", n)
    incomplete <- c("Q1", "Q4", "Q6", "Q11", "Q12", "D1", "D2", "D5")
    for (j in seq_along(incomplete)) d[[incomplete[j]]][c(j, j + 80L)] <- NA
    d
  })
  original <- mie_encode(raw)
  mie_validate(original)
  stopifnot(nrow(original) == nrow(raw), identical(original$Observation, raw$Observation),
            all(original$D1[!is.na(original$D1)] %in% 1:5),
            all(original$D5[!is.na(original$D5)] %in% 1:5),
            all(original$D2[!is.na(original$D2)] %in% 0:1),
            identical(is.na(original$D1), is.na(raw$D1)),
            identical(is.na(original$D2), is.na(raw$D2)),
            identical(is.na(original$D5), is.na(raw$D5)))
  scenarios <- mie_scenarios()
  stopifnot(nrow(scenarios) == 5L, !anyDuplicated(scenarios$scenario),
            sum(scenarios$delta_B == 0 & scenarios$delta_A == 0) == 1L,
            setequal(paste(sign(scenarios$delta_B), sign(scenarios$delta_A)),
                     c("0 0", "-1 -1", "1 1", "-1 1", "1 -1")),
            all(abs(scenarios$delta_B) %in% c(0, log(2))),
            all(abs(scenarios$delta_A) %in% c(0, log(2))))
  baseline <- scenarios[scenarios$delta_B == 0 & scenarios$delta_A == 0, , drop = FALSE]
  initial <- mie_initialize(original, streams[[5L]], baseline, 20261111L, 1L)
  stopifnot(!anyNA(initial$data[, expected_variables, drop = FALSE]))
  for (target in c("Q6", "D2")) {
    # Start with n rows, then form indicator columns independently of the
    # production model-matrix code and its chosen column labels.
    expected_design <- matrix(1, nrow(original), 1L)
    for (item in setdiff(expected_variables, target)) {
      indicators <- if (item == "D2") 1L else 2:5
      for (level in indicators)
        expected_design <- cbind(expected_design, as.numeric(initial$data[[item]] == level))
    }
    for (site in 2:8)
      expected_design <- cbind(expected_design, as.numeric(as.character(raw$Site) == site))
    expected_design <- cbind(expected_design, as.numeric(raw$Language == "Spanish"))
    actual_design <- mie_design(initial$data, target)
    close(unname(actual_design), unname(expected_design))
    stopifnot(ncol(actual_design) == if (target == "D2") 93L else 90L,
              !anyDuplicated(colnames(actual_design)),
              !any(grepl("D8", colnames(actual_design))),
              sum(colSums(abs(actual_design)) == 0) >= 4L)
  }
  passed("auxiliary encoding, all five scenarios, fixed reference coding, and retained zero columns")

  # Reconstruct the matching operation on the same stream. This checks the
  # distinction between donor posterior means and recipient posterior draws,
  # target scaling, binary support, and use of the staffing tilt only when due.
  scenario <- scenarios[scenarios$scenario == "busy_lower", , drop = FALSE]
  for (target in c("Q1", "Q11", "D2")) {
    observed <- which(!is.na(original[[target]]))
    missing <- which(is.na(original[[target]]))
    target_X <- mie_design(initial$data, target)
    scaled <- original[[target]][observed]
    if (target != "D2") scaled <- (scaled - 1) / 4
    delta <- if (target == "Q11") scenario$delta_B else 0
    expected_selected <- mie_with_rng(streams[[6L]], {
      posterior <- mie_posterior(target_X[observed, , drop = FALSE], scaled)
      match <- mie_draw_donors(
        as.numeric(target_X[observed, , drop = FALSE] %*% posterior$mean),
        as.numeric(target_X[missing, , drop = FALSE] %*% posterior$draw),
        original[[target]][observed], delta)
      observed[match$selected]
    })
    updated <- mie_with_rng(streams[[6L]], mie_update(initial$data, original, target, scenario))
    stopifnot(identical(updated$selected, expected_selected),
              identical(updated$data[[target]][missing], original[[target]][expected_selected]),
              identical(updated$data[[target]][observed], original[[target]][observed]))
  }
  passed("Type-1 matching predictions, ordinal/binary scaling, and original-value donor copying")

  constant_target <- original
  constant_target$Q27[1L] <- NA_integer_
  mie_validate(constant_target)
  constant_initial <- mie_initialize(constant_target, streams[[6L]], baseline, 20261111L, 1L)
  constant_update <- mie_with_rng(streams[[6L]],
    mie_update(constant_initial$data, constant_target, "Q27", baseline))
  stopifnot(all(constant_update$data$Q27 == 1L),
            !is.na(constant_target$Q27[constant_update$selected]),
            all(constant_update$details$unique_categories == 1L),
            all(abs(constant_update$details$conditional_mean_shift) < 1e-12))
  passed("a constant observed target remains valid under the proper prior and donor support")

  # Checkpoint resumption must preserve the exact random stream and imputations.
  direct <- mie_advance(initial, original, 3L)
  resumed <- mie_advance(mie_advance(initial, original, 1L), original, 3L)
  repeated <- mie_advance(mie_initialize(original, streams[[5L]], baseline,
                                        20261111L, 1L), original, 3L)
  stopifnot(direct$status == "complete", direct$iteration == 3L,
            identical(direct$data, resumed$data), identical(direct$stream, resumed$stream),
            identical(direct$data, repeated$data), identical(direct$stream, repeated$stream),
            identical(direct$trace, resumed$trace), identical(direct$donor_counts, resumed$donor_counts),
            identical(direct$staffing_updates, resumed$staffing_updates))
  for (item in expected_variables) {
    observed <- !is.na(original[[item]])
    support <- if (item == "D2") 0:1 else 1:5
    stopifnot(identical(direct$data[[item]][observed], original[[item]][observed]),
              all(direct$data[[item]] %in% support),
              all(direct$data[[item]][!observed] %in% original[[item]][observed]))
    if (any(!observed)) {
      stopifnot(sum(direct$donor_counts[[item]]) == 3L * sum(!observed),
                all(direct$donor_counts[[item]][!observed] == 0L))
      mapping <- direct$terminal_donors[[item]]
      donors <- match(mapping$selected_record, original$Observation)
      recipients <- match(mapping$record, original$Observation)
      stopifnot(!anyNA(donors), !anyNA(recipients), all(observed[donors]),
                all(!observed[recipients]),
                all(direct$data[[item]][recipients] == original[[item]][donors]))
    }
  }
  stopifnot(identical(direct$data$Observation, original$Observation),
            identical(direct$data$Site, original$Site),
            identical(direct$data$Language, original$Language))
  for (outcome in c(6L, 4L)) {
    analysis <- mie_analysis_data(original, direct$data, outcome)
    observed <- !is.na(original[[paste0("Q", outcome)]])
    stopifnot(identical(analysis$Observation, original$Observation[observed]),
              !anyNA(analysis[[paste0("Q", outcome)]]),
              nrow(analysis) > sum(complete.cases(original[, expected_variables, drop = FALSE])))
  }
  passed("exact resumption, observed values and donor eligibility, support, and outcome-observed targets")

  bad_site <- original
  bad_site$Site[1L] <- NA
  expect_error(mie_validate(bad_site))
  bad_language <- original
  bad_language$Language[1L] <- NA
  expect_error(mie_validate(bad_language))
  bad_value <- original
  bad_value$Q11[1L] <- 6L
  expect_error(mie_validate(bad_value))
  bad_id <- original
  bad_id$Observation[2L] <- bad_id$Observation[1L]
  expect_error(mie_validate(bad_id))
  passed("missing fixed predictors, invalid ordinal support, and duplicate identifiers fail validation")

  damaged <- initial
  damaged$data$Q2[1L] <- NA_integer_
  failed_chain <- mie_advance(damaged, original, 1L)
  stopifnot(failed_chain$status == "failed", length(failed_chain$failure$reason) == 1L,
            identical(failed_chain, mie_advance(failed_chain, original, 3L)))
  passed("conditional failures are retained and cannot be silently resumed or replaced")

  # Summary tests use fixed coefficient lists, never an ordinal outcome fit.
  registry <- list()
  for (scenario in scenarios$scenario) for (seed in c(20261111L, 20261112L))
    for (outcome in c(6L, 4L)) for (chain in 1:2) {
      values <- c(beta_B = chain * 0.02, beta_A = chain * 0.03)
      if (outcome == 6L) values <- c(values, d = unname(values["beta_A"] - values["beta_B"]))
      registry[[length(registry) + 1L]] <- list(
        scenario = scenario, seed = seed, chain = as.integer(chain), iteration = 3L,
        outcome = outcome, n = 198L, valid = TRUE, reason = "", values = values)
    }
  summary <- mie_summarize_fits(registry, m = 2L)
  stopifnot(nrow(summary$estimates) == 50L,
            all(summary$availability$status == "available"),
            all(summary$availability$attempted == 2L),
            !any(grepl("p_value|pvalue|standard_error|ci_lower|ci_upper", names(summary$estimates))))
  difference <- summary$estimates[summary$estimates$component == "d", , drop = FALSE]
  close(difference$estimate, rep(0.015, nrow(difference)))
  close(difference$exponentiated, rep(exp(0.015), nrow(difference)))
  stopifnot(all(difference$exponentiated_quantity == "ratio_of_odds_ratios"))

  damaged_registry <- registry
  failed_index <- which(vapply(damaged_registry, function(x)
    x$scenario == "baseline" && x$seed == 20261111L && x$outcome == 6L && x$chain == 1L,
    logical(1L)))
  damaged_registry[[failed_index]]$valid <- FALSE
  damaged_registry[[failed_index]]$reason <- "synthetic retained failure"
  damaged_registry[[failed_index]]$values[] <- NA_real_
  failed_summary <- mie_summarize_fits(damaged_registry, m = 2L)
  stopifnot(!any(failed_summary$estimates$scenario == "baseline" & failed_summary$estimates$outcome == 6L),
            any(failed_summary$estimates$scenario == "baseline" & failed_summary$estimates$outcome == 4L),
            sum(failed_summary$availability$failed) == 1L)
  missing_summary <- mie_summarize_fits(registry[-failed_index], m = 2L)
  stopifnot(!any(missing_summary$estimates$scenario == "baseline" & missing_summary$estimates$outcome == 6L))
  duplicate_registry <- registry
  duplicate_registry[[failed_index + 1L]]$chain <- 1L
  duplicate_summary <- mie_summarize_fits(duplicate_registry, m = 2L)
  stopifnot(!any(duplicate_summary$estimates$scenario == "baseline" & duplicate_summary$estimates$outcome == 6L))
  passed("log-scale averaging, Q6 contrast, and suppression after failed, missing, or duplicate fits")

  # The remaining checks create synthetic files only inside this test-owned
  # temporary directory. They exercise scheduling, provenance, and review
  # boundaries without fitting an ordinal outcome model.
  test_root <- tempfile("mie-extension-tests-", tmpdir = tempdir())
  dir.create(test_root, recursive = TRUE)
  on.exit({
    resolved <- normalizePath(test_root, winslash = "/", mustWork = FALSE)
    temporary_parent <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
    if (!startsWith(resolved, paste0(temporary_parent, "/mie-extension-tests-")))
      stop("Refusing cleanup outside the test-owned temporary directory.")
    unlink(resolved, recursive = TRUE)
  }, add = TRUE)
  one_directory <- file.path(test_root, "one_worker")
  two_directory <- file.path(test_root, "two_workers")
  dir.create(one_directory)
  dir.create(two_directory)
  parallel_streams <- setNames(lapply(c(20261111L, 20261112L), function(seed)
    mie_streams(seed, 2L)), as.character(c(20261111L, 20261112L)))
  tasks_one <- mie_tasks(one_directory, m = 2L, iterations = 2L, scenarios = "baseline")
  tasks_two <- mie_tasks(two_directory, m = 2L, iterations = 2L, scenarios = "baseline")
  one_results <- mie_parallel(tasks_one, mie_generation_worker, script_dir, original,
                              parallel_streams, workers = 1L)
  two_results <- mie_parallel(tasks_two, mie_generation_worker, script_dir, original,
                              parallel_streams, workers = 2L)
  stopifnot(all(vapply(one_results, function(x) identical(x$status, "complete"), logical(1L))),
            all(vapply(two_results, function(x) identical(x$status, "complete"), logical(1L))))
  for (i in seq_len(nrow(tasks_one))) {
    first <- readRDS(tasks_one$path[i])
    second <- readRDS(tasks_two$path[i])
    stopifnot(identical(first, second))
  }
  passed("one-worker and two-worker PSOCK runs preserve identical chain states and streams")

  diagnostics <- mie_write_diagnostics(one_directory, original, m = 2L,
                                       iterations = 2L, scenarios = "baseline")
  stopifnot(isTRUE(diagnostics$complete))
  mie_validate_diagnostics(diagnostics)
  regenerated <- mie_write_diagnostics(one_directory, original, m = 2L,
                                       iterations = 2L, scenarios = "baseline")
  stopifnot(identical(diagnostics$fingerprint, regenerated$fingerprint))
  mie_validate_diagnostics(regenerated)
  passed("diagnostic regeneration preserves its fingerprint and content checks")

  manifest_path <- file.path(one_directory, "synthetic_manifest.rds")
  review_path <- file.path(one_directory, "synthetic_review.rds")
  saveRDS(list(version = "1.4", scenario = "baseline", m = 2L, iterations = 2L), manifest_path)
  expect_error(mie_require_review(manifest_path, diagnostics, review_path))
  review <- list(decision = "accept", reviewer = "synthetic fixture", note = "Synthetic diagnostics only.",
                 manifest_md5 = unname(tools::md5sum(manifest_path)),
                 diagnostic_fingerprint = diagnostics$fingerprint)
  saveRDS(review, review_path)
  accepted <- mie_require_review(manifest_path, diagnostics, review_path)
  stopifnot(identical(accepted$decision, "accept"))
  rejected <- review
  rejected$decision <- "extend_iterations"
  saveRDS(rejected, review_path)
  expect_error(mie_require_review(manifest_path, diagnostics, review_path))
  retained <- mie_require_review(manifest_path, diagnostics, review_path, require_accept = FALSE)
  stopifnot(identical(retained$decision, "extend_iterations"))
  saveRDS(review, review_path)
  incomplete_diagnostics <- diagnostics
  incomplete_diagnostics$complete <- FALSE
  expect_error(mie_require_review(manifest_path, incomplete_diagnostics, review_path))

  # Editing a diagnostic after review must invalidate the fit gate, even if the
  # review file and top-level index have not changed. Restore original bytes.
  trace_path <- file.path(diagnostics$directory, "trace.csv")
  original_bytes <- readBin(trace_path, what = "raw", n = file.info(trace_path)$size)
  cat("tampered diagnostic\n", file = trace_path, append = TRUE)
  expect_error(mie_require_review(manifest_path, diagnostics, review_path))
  writeBin(original_bytes, trace_path)
  mie_require_review(manifest_path, diagnostics, review_path)
  changed_review <- review
  changed_review$manifest_md5 <- "not_the_manifest_digest"
  saveRDS(changed_review, review_path)
  expect_error(mie_require_review(manifest_path, diagnostics, review_path))
  saveRDS(review, review_path)
  passed("generation diagnostics require a current accepted review before the fit gate opens")

  fit_tasks <- mie_fit_tasks(one_directory, m = 2L, iterations = 2L, scenarios = "baseline")
  fit_task <- as.list(fit_tasks[1L, , drop = FALSE])
  stopifnot(is.null(mie_load_saved_fit(fit_task)))
  cached <- list(scenario = fit_task$scenario, seed = fit_task$seed, chain = fit_task$chain,
                 iteration = fit_task$iterations, outcome = fit_task$outcome, n = 198L,
                 valid = TRUE, reason = "", values = c(beta_B = 0.2, beta_A = 0.3, d = 0.1))
  dir.create(dirname(fit_task$fit_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cached, fit_task$fit_path)
  stopifnot(identical(mie_load_saved_fit(fit_task), cached))
  retained_failure <- cached
  retained_failure$valid <- FALSE
  retained_failure$reason <- "synthetic retained fit failure"
  retained_failure$values[] <- NA_real_
  sidecar <- paste0(fit_task$fit_path, ".failure.rds")
  saveRDS(retained_failure, sidecar)
  cached_digest <- unname(tools::md5sum(fit_task$fit_path))
  stopifnot(identical(mie_load_saved_fit(fit_task), retained_failure),
            identical(unname(tools::md5sum(fit_task$fit_path)), cached_digest))
  # The retained failure must still win when the cached main file is unreadable.
  writeBin(charToRaw("synthetic corrupt cache"), fit_task$fit_path)
  stopifnot(identical(mie_load_saved_fit(fit_task), retained_failure))
  wrong_identity <- retained_failure
  wrong_identity$chain <- wrong_identity$chain + 1L
  saveRDS(wrong_identity, sidecar)
  expect_error(mie_load_saved_fit(fit_task))
  saveRDS(retained_failure, sidecar)
  passed("summary reload prioritizes retained failures, preserves caches, and rejects identity mismatches")

  data.frame(test = checks, passed = TRUE, detail = "Synthetic check passed.", stringsAsFactors = FALSE)
}
