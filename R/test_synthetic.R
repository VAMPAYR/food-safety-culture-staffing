# Synthetic verification only. Sourcing this file does not run tests.
# After independent pre-execution review, call:
# run_synthetic_tests(".../R/model_engine.R", output_dir = ".../verification")

run_synthetic_tests <- function(engine_path, output_dir = NULL, seed = 61723L) {
  source(engine_path, local = environment())
  set.seed(seed)
  results <- list()
  check <- function(name, pass, detail = "") {
    results[[length(results) + 1L]] <<- data.frame(test = name,
      passed = isTRUE(pass), detail = as.character(detail), stringsAsFactors = FALSE)
  }
  n <- 2400L
  d <- data.frame(Observation = seq_len(n), Site = rep(1:8, each = n / 8))
  for (v in c("L_c", "E_c", "R_c", "M_c", "B_c", "A_c")) d[[v]] <- rnorm(n, sd = 0.7)
  site_effect <- c(0, -0.2, 0.1, 0.15, -0.1, 0.2, -0.15, 0.05)
  eta <- site_effect[d$Site] + 0.25 * d$L_c + 0.15 * d$E_c +
    0.2 * d$R_c - 0.15 * d$M_c + 0.45 * d$B_c + 0.7 * d$A_c
  cut <- c(-1.2, -0.15, 1.1)
  F <- vapply(cut, function(t) plogis(t - eta), numeric(n))
  u <- runif(n)
  y <- 1L + rowSums(u > F)
  d$Y <- ordered(y, levels = 1:4)
  predictors <- c("L_c", "E_c", "R_c", "M_c", "B_c", "A_c")
  fit <- fit_ordinal(d, "Y", predictors, model_id = "synthetic_full")
  check("ordinal_fit_accepted", fit$valid, fit$status)
  check("eight_sites_reference_available", fit$inference_valid)
  if (!length(fit$coef)) {
    out <- do.call(rbind, results)
    if (!is.null(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      write.csv(out, file.path(output_dir, "synthetic_tests.csv"), row.names = FALSE)
      saveRDS(fit, file.path(output_dir, "synthetic_failed_fit.rds"))
    }
    stop("Synthetic ordinal fit returned no parameters; inspect preserved results.")
  }
  check("positive_sign_convention", fit$coef["B_c"] > 0 && fit$coef["A_c"] > 0)
  check("known_positive_coefficients_recovered", abs(fit$coef["B_c"] - 0.45) < 0.2 &&
      abs(fit$coef["A_c"] - 0.7) < 0.2)
  check("category_probability_sums", max(abs(rowSums(fit$probabilities) - 1)) < 1e-8)
  check("score_tolerance", fit$diagnostics$max_score_per_n < 1e-6)
  check("all_parameters_included", ncol(fit$scores) == 16L && nrow(fit$H) == 16L &&
      sum(grepl("^threshold:", names(fit$coef))) == 3L)

  # Independently evaluated probabilities and worker log likelihoods. These
  # use direct differences on moderate synthetic parameters, without the
  # engine's stable-probability or score routines.
  X <- fit$design; p <- ncol(X); theta <- fit$coef; q <- length(theta)
  independent_loglik <- function(b, individual = FALSE) {
    et <- as.vector(X %*% b[seq_len(p)])
    tt <- b[p + seq_len(3L)]
    ff <- vapply(tt, function(t) 1 / (1 + exp(-(t - et))), numeric(n))
    pp <- cbind(ff[, 1L], ff[, 2L] - ff[, 1L], ff[, 3L] - ff[, 2L], 1 - ff[, 3L])
    lli <- log(pp[cbind(seq_len(n), y)])
    if (individual) lli else sum(lli)
  }
  step_score <- 1e-5
  numerical_scores <- matrix(0, n, q)
  for (j in seq_len(q)) {
    shift <- numeric(q); shift[j] <- step_score
    numerical_scores[, j] <- (independent_loglik(theta + shift, TRUE) -
      independent_loglik(theta - shift, TRUE)) / (2 * step_score)
  }
  check("analytic_worker_scores_match_finite_differences",
    max(abs(numerical_scores - fit$scores)) < 2e-7,
    max(abs(numerical_scores - fit$scores)))
  step_hessian <- 2e-4
  H_numeric <- matrix(0, q, q); center_ll <- independent_loglik(theta)
  for (j in seq_len(q)) for (k in seq_len(j)) {
    ej <- ek <- numeric(q); ej[j] <- step_hessian; ek[k] <- step_hessian
    if (j == k) {
      value <- -(independent_loglik(theta + ej) - 2 * center_ll +
        independent_loglik(theta - ej)) / step_hessian^2
    } else {
      value <- -(independent_loglik(theta + ej + ek) - independent_loglik(theta + ej - ek) -
        independent_loglik(theta - ej + ek) + independent_loglik(theta - ej - ek)) /
        (4 * step_hessian^2)
    }
    H_numeric[j, k] <- H_numeric[k, j] <- value
  }
  check("observed_information_matches_numerical_hessian", max(abs(H_numeric - fit$H)) < 0.003,
    max(abs(H_numeric - fit$H)))
  check("likelihood_matches_independent_probabilities", abs(center_ll - fit$logLik) < 1e-8)
  # Independently aggregate numeric worker scores using explicit site loops.
  S_numeric <- do.call(rbind, lapply(1:8, function(g) colSums(numerical_scores[d$Site == g, , drop = FALSE])))
  check("site_aggregation_matches_independent_scores", max(abs(S_numeric - fit$site_scores)) < 1e-5)
  bread_numeric <- solve(H_numeric)
  meat_numeric <- Reduce(`+`, lapply(seq_len(8L), function(g) tcrossprod(S_numeric[g, ])))
  V_numeric <- (8 / 7) * bread_numeric %*% meat_numeric %*% bread_numeric
  check("site_sandwich_matches_independent_matrix_calculation", max(abs(V_numeric - fit$V)) < 2e-6,
    max(abs(V_numeric - fit$V)))
  # At the accepted finite numerical optimum the remaining eighth singular
  # value can reflect score tolerance; assess it relative to the leading one.
  sv <- svd(fit$V)$d
  check("site_sandwich_rank_limit", sv[8L] / sv[1L] < 1e-7,
    paste("eighth/first singular value", sv[8L] / sv[1L]))
  check("full_sandwich_singular", sv[length(sv)] / sv[1L] < 1e-12)

  difference <- contrast_summary(fit, c(A_c = 1, B_c = -1), "A_minus_B")
  variance <- fit$V["A_c", "A_c"] + fit$V["B_c", "B_c"] - 2 * fit$V["A_c", "B_c"]
  check("difference_contrast_covariance", difference$valid && abs(difference$se^2 - variance) < 1e-12)
  check("difference_uses_t7_reference", abs(difference$p_reference -
      2 * pt(-abs(difference$estimate / sqrt(variance)), 7)) < 1e-12)
  check("difference_interval_uses_t7", abs(difference$upper - difference$estimate -
      qt(0.975, 7) * sqrt(variance)) < 1e-12)
  joint <- joint_summary(fit, c("B_c", "A_c"), "staffing_joint")
  bv <- fit$coef[c("B_c", "A_c")]; vv <- fit$V[c("B_c", "A_c"), c("B_c", "A_c")]
  W <- as.numeric(t(bv) %*% solve(vv) %*% bv)
  check("joint_uses_F2_7_reference", joint$valid && abs(joint$statistic - W / 2) < 1e-10 &&
      abs(joint$p_reference - pf(W / 2, 2, 7, lower.tail = FALSE)) < 1e-12)
  bad_cov <- fit; bad_cov$V["A_c", ] <- bad_cov$V["B_c", ]
  bad_cov$V[, "A_c"] <- bad_cov$V[, "B_c"]
  bad_joint <- joint_summary(bad_cov, c("B_c", "A_c"))
  check("singular_contrast_rejected", !bad_joint$valid && is.na(bad_joint$p_reference))
  check("failed_secondary_member_blocks_entire_Holm_family",
    all(is.na(secondary_holm(c(RQ2 = 0.01, RQ3 = NA_real_)))))
  check("Holm_matches_prespecified_formula", all.equal(unname(secondary_holm(c(0.04, 0.01))),
    c(0.04, 0.02)) == TRUE)

  check("alternate_fit_staffing_agreement", max(abs(fit$coef[c("B_c", "A_c")] -
    fit$alternate[c("B_c", "A_c")])) <= 1e-4)
  alt_ll <- independent_loglik(fit$alternate)
  check("alternate_fit_likelihood_agreement", abs(alt_ll - fit$logLik) <= 1e-6 * (1 + abs(fit$logLik)))

  # Deliberately stop a real synthetic polr fit at a loose objective tolerance.
  # The numerical refinement must preserve that failed score attempt and
  # converge under the unchanged acceptance criterion on the same model.
  retry_prep <- .prepare_model(d, "Y", predictors, "Site")
  refined <- NULL
  loose <- .capture(MASS::polr(retry_prep$formula, data = retry_prep$data,
    method = "logistic", Hess = FALSE, control = list(maxit = 5000L, reltol = 1e-4)))
  if (is.null(loose$value)) {
    check("score_refinement_fixture_fit_available", FALSE, loose$error)
  } else {
    loose_theta <- c(coef(loose$value)[colnames(retry_prep$X)],
      setNames(loose$value$zeta, paste0("threshold:", names(loose$value$zeta))))
    original_attempt <- list(primary = list(solver = "MASS::polr loose synthetic fixture",
      warnings = loose$warnings, error = loose$error, convergence = loose$value$convergence,
      theta = loose_theta))
    refined <- .refine_ordinal_primary(retry_prep, loose_theta, loose$value$convergence, original_attempt)
    check("score_refinement_fixture_exceeds_acceptance", loose$value$convergence == 0 &&
      refined$attempts$primary$max_score_per_n >= 1e-6)
    check("score_refinement_preserves_initial_attempt", identical(refined$attempts$primary$theta, loose_theta) &&
      identical(refined$attempts$primary$convergence, loose$value$convergence) &&
      identical(refined$attempts$primary$warnings, loose$warnings))
    check("score_refinement_converges_without_relaxing_acceptance",
      identical(refined$attempts$selected_primary, "primary_refinement") &&
      refined$convergence == 0 && max(abs(refined$comp$score)) / n < 1e-6)
    check("score_refinement_agrees_with_independent_accepted_fit",
      max(abs(refined$theta[c("B_c", "A_c")] - fit$alternate[c("B_c", "A_c")])) <= 1e-4 &&
      abs(refined$comp$logLik - fit$logLik) <= 1e-6 * (1 + abs(fit$logLik)))
    check("score_refinement_preserves_likelihood", isTRUE(refined$attempts$primary_refinement$logLik_change >=
      -refined$attempts$primary_refinement$allowed_logLik_decrease))
    no_retry <- .refine_ordinal_primary(retry_prep, fit$coef, 0L,
      list(primary = list(theta = fit$coef, convergence = 0L)))
    check("accepted_score_does_not_trigger_refinement", is.null(no_retry$attempts$primary_refinement))
  }

  reverse <- d; reverse$M_c <- -reverse$M_c
  rev_fit <- fit_ordinal(reverse, "Y", predictors, model_id = "synthetic_global_M_reversal")
  check("global_management_reversal_fit_valid", rev_fit$valid, rev_fit$status)
  check("global_management_reversal_staffing_invariance", max(abs(fit$coef[c("B_c", "A_c")] -
      rev_fit$coef[c("B_c", "A_c")])) < 1e-4)
  check("global_management_reversal_slope_sign", abs(fit$coef["M_c"] + rev_fit$coef["M_c"]) < 1e-4)
  check("global_management_reversal_probability_invariance",
    max(abs(fit$probabilities - rev_fit$probabilities)) < 1e-4)
  check("global_management_reversal_likelihood_invariance", abs(fit$logLik - rev_fit$logLik) < 1e-5)

  reduced <- fit_ordinal(d, "Y", setdiff(predictors, c("B_c", "A_c")), model_id = "synthetic_reduced")
  comparison <- same_row_comparison(fit, reduced)
  check("same_row_nested_comparison", comparison$valid && comparison$delta_logLik >= 0)
  wrong_rows <- reduced; wrong_rows$rows <- rev(wrong_rows$rows)
  check("different_rows_comparison_rejected", !same_row_comparison(fit, wrong_rows)$valid)
  rank_bad <- d; rank_bad$duplicate <- rank_bad$B_c
  rejected <- fit_ordinal(rank_bad, "Y", c(predictors, "duplicate"), model_id = "synthetic_rank_deficient")
  check("rank_deficient_design_rejected", !rejected$valid && grepl("rank", rejected$status))

  binary_fits <- lapply(1:3, function(k) {
    dd <- d; dd$Z <- as.integer(y > k)
    fit_binary(dd, "Z", predictors, model_id = paste0("synthetic_threshold_", k))
  })
  check("binary_threshold_diagnostics_valid", all(vapply(binary_fits, function(f) f$valid, logical(1))))
  check("binary_higher_response_sign", all(vapply(binary_fits, function(f)
    f$coef["B_c"] > 0 && f$coef["A_c"] > 0, logical(1))))
  bf <- binary_fits[[2L]]; by <- as.integer(y > 2L)
  numeric_binary_score <- bf$design * as.vector(by - plogis(drop(bf$design %*% bf$coef)))
  check("binary_scores_match_direct_formula", max(abs(numeric_binary_score - bf$scores)) < 1e-12)
  deleted <- fit_ordinal(d[d$Site != 1L, ], "Y", predictors, model_id = "synthetic_site_1_deleted")
  check("site_deletion_keeps_point_fit_only", deleted$valid && !deleted$inference_valid && deleted$G == 7L,
    deleted$status)
  deleted_summary <- contrast_summary(deleted, c(B_c = 1))
  check("site_deletion_has_no_reference_interval", is.finite(deleted_summary$estimate) &&
      is.na(deleted_summary$lower) && is.na(deleted_summary$p_reference))

  categorical <- d
  for (v in c("B", "A")) {
    x <- d[[paste0(v, "_c")]]
    categorical[[paste0(v, "_category")]] <- factor(ifelse(x < -0.3, "disagreement",
      ifelse(x > 0.3, "agreement", "neutrality")),
      levels = c("neutrality", "disagreement", "agreement"))
  }
  category_fit <- fit_ordinal(categorical, "Y",
    c(setdiff(predictors, c("B_c", "A_c")), "B_category", "A_category"),
    model_id = "synthetic_categorical_staffing")
  category_names <- c("B_categorydisagreement", "B_categoryagreement",
    "A_categorydisagreement", "A_categoryagreement")
  check("categorical_fit_and_parameter_labels", category_fit$valid &&
      all(category_names %in% names(category_fit$coef)), category_fit$status)
  check("categorical_staffing_alternate_agreement", all(category_names %in% names(category_fit$alternate)) &&
      max(abs(category_fit$coef[category_names] - category_fit$alternate[category_names])) <= 1e-4)

  # Arithmetic unit fixtures supply accepted deletion-fit records with known
  # staffing perturbations. They test the covariance wrapper, not fitting.
  perturbation <- cbind(B_c = c(-3, -2, -1, 0, 1, 2, 3, 4) / 100,
    A_c = c(1, -2, 3, -4, 2, -1, 4, -3) / 100)
  make_deletion_fixtures <- function(base) lapply(1:8, function(g) {
    f <- base; keep <- as.character(base$site) != as.character(g)
    f$site <- droplevels(base$site[keep]); f$rows <- base$rows[keep]
    f$outcome <- base$outcome[keep]; f$n <- sum(keep); f$G <- 7L
    f$coef[c("B_c", "A_c")] <- base$coef[c("B_c", "A_c")] + perturbation[g, ]
    f$inference_valid <- FALSE; f
  })
  deletion_fixtures <- make_deletion_fixtures(fit)
  jk <- cluster_jackknife_staffing(fit, deletion_fixtures, "Q6")
  jk_hand <- matrix(0, 2L, 2L)
  for (g in 1:8) for (j in 1:2) for (k in 1:2)
    jk_hand[j, k] <- jk_hand[j, k] + (7 / 8) * perturbation[g, j] * perturbation[g, k]
  check("jackknife_covariance_matches_hand_calculation", all(jk$valid) &&
      max(abs(attr(jk, "Vjk") - jk_hand)) < 1e-14)
  mean_centered <- (7 / 8) * crossprod(sweep(perturbation, 2, colMeans(perturbation)))
  check("jackknife_uses_full_fit_center", max(abs(attr(jk, "Vjk") - mean_centered)) > 1e-6)
  jk_d <- jk[jk$term == "A_minus_B", ]
  jk_d_variance <- jk_hand[1, 1] + jk_hand[2, 2] - 2 * jk_hand[1, 2]
  check("jackknife_difference_includes_covariance", abs(jk_d$se^2 - jk_d_variance) < 1e-14)
  check("jackknife_interval_uses_t7", abs(jk_d$upper - jk_d$estimate -
      qt(0.975, 7) * sqrt(jk_d_variance)) < 1e-12)
  check("jackknife_omits_reference_tests", !any(grepl("p_reference|statistic|Holm", names(jk))))
  jk_q4 <- cluster_jackknife_staffing(fit, deletion_fixtures, "Q4")
  check("jackknife_Q4_has_no_difference", identical(jk_q4$term, c("B_c", "A_c")))
  reversed_jk <- cluster_jackknife_staffing(rev_fit, make_deletion_fixtures(rev_fit), "Q6")
  check("jackknife_global_management_reversal_invariance", all(reversed_jk$valid) &&
      max(abs(attr(jk, "Vjk") - attr(reversed_jk, "Vjk"))) < 1e-12 &&
      max(abs(jk$estimate - reversed_jk$estimate)) < 1e-4)
  reordered_jk <- cluster_jackknife_staffing(fit, rev(deletion_fixtures), "Q6")
  check("jackknife_deletion_order_invariance", max(abs(attr(jk, "Vjk") - attr(reordered_jk, "Vjk"))) < 1e-14)
  failed_deletions <- deletion_fixtures; failed_deletions[[3L]]$valid <- FALSE
  failed_deletions[[3L]]$status <- "synthetic_failure"
  failed_jk <- cluster_jackknife_staffing(fit, failed_deletions, "Q6")
  check("jackknife_failed_deletion_blocks_entire_outcome", !any(failed_jk$valid) &&
      all(is.na(failed_jk$se)) && all(failed_jk$n_valid_deletions == 7L) &&
      "synthetic_failure" %in% attr(failed_jk, "deletion_status"))
  duplicated_deletions <- deletion_fixtures; duplicated_deletions[[8L]] <- duplicated_deletions[[7L]]
  duplicate_jk <- cluster_jackknife_staffing(fit, duplicated_deletions, "Q6")
  check("jackknife_duplicate_deleted_site_rejected", !any(duplicate_jk$valid))

  # Check stable likelihood and score evaluations in tails and narrow gaps,
  # independent of fitting or accepting an extreme model.
  X_extreme <- matrix(c(-1000, 1000, 0), ncol = 1L)
  extreme <- ordinal_components(c(slope = 1, t1 = -1, t2 = 0, t3 = 1), X_extreme,
    c(4L, 1L, 2L), 4L)
  check("extreme_tail_loglikelihood_finite", is.finite(extreme$logLik) && all(is.finite(extreme$scores)))
  narrow <- ordinal_components(c(slope = 0, t1 = -1, t2 = -1 + 1e-8, t3 = 1),
    matrix(1, 1L, 1L), 2L, 4L)
  check("narrow_threshold_loglikelihood_finite", is.finite(narrow$logLik) && all(is.finite(narrow$scores)))
  out <- do.call(rbind, results)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(out, file.path(output_dir, "synthetic_tests.csv"), row.names = FALSE)
    saveRDS(list(results = out, seed = seed, full = fit, reversed = rev_fit,
      reduced = reduced, binary = binary_fits, deleted = deleted,
      categorical = category_fit, jackknife = jk, jackknife_failed = failed_jk,
      score_refinement = refined,
      independent = list(scores = numerical_scores, H = H_numeric, V = V_numeric),
      session = sessionInfo()), file.path(output_dir, "synthetic_verification.rds"))
  }
  if (any(!out$passed)) stop("Synthetic verification failed: ", paste(out$test[!out$passed], collapse = ", "))
  out
}
