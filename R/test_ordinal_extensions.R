# Synthetic verification for specification 1.4 priority 2.
# Sourcing defines functions only. Do not execute before the independent gate.
# run_ordinal_extension_tests(script_dir, output_dir = NULL, seed = 20261114L)

run_ordinal_extension_tests <- function(script_dir, output_dir = NULL, seed = 20261114L) {
  for (nm in c("data_prep.R", "model_engine.R", "ordinal_extensions.R"))
    source(file.path(script_dir, nm), local = environment())
  if (!requireNamespace("ordinal", quietly = TRUE)) stop("ordinal is needed for verification.")
  set.seed(seed)
  records <- list(); saved <- list()
  check <- function(name, pass, detail = "") {
    records[[length(records) + 1L]] <<- data.frame(test = name, passed = isTRUE(pass),
      detail = as.character(detail), stringsAsFactors = FALSE)
  }
  scaled_error <- function(a, b) max(abs(a - b)) / (1 + max(abs(b)))
  n <- 1600L
  d <- data.frame(Observation = as.character(seq_len(n)), Site = rep(as.character(1:8), each = n / 8L))
  for (v in c("L_c", "E_c", "R_c", "M_c")) d[[v]] <- runif(n, -0.9, 0.9)
  d$B <- sample(1:5, n, TRUE); d$A <- sample(1:5, n, TRUE)
  d$B_c <- d$B - 3; d$A_c <- d$A - 3
  d$Y <- ordered(rep(1:4, length.out = n), levels = 1:4)
  d <- .oe_data(d)
  z <- oe_ppo_design(d)
  theta <- setNames(numeric(20L), z$names)
  theta[1:3] <- c(-1.1, 0.3, 1.7)
  theta[4:6] <- c(0.16, 0.22, 0.28); theta[7:9] <- c(0.25, 0.21, 0.17)
  theta[c("L_c", "M_c", "B_c", "A_c")] <- c(0.2, -0.1, 0.25, 0.4)
  theta[grep("^Site", names(theta))] <- seq(-0.15, 0.15, length.out = 7)
  c0 <- oe_ppo_components(theta, z)
  u <- runif(n)
  d$Y <- ordered(1L + rowSums(u > t(apply(c0$probabilities, 1, cumsum))[, 1:3]), levels = 1:4)
  z <- oe_ppo_design(d)
  comp <- oe_ppo_components(theta, z)
  check("PPO_dimensions", z$q == 20L && ncol(comp$scores) == 20L && nrow(comp$H) == 20L)
  check("PPO_probabilities_sum_to_one", max(abs(rowSums(comp$probabilities) - 1)) < 1e-12)
  check("PPO_strict_probabilities", all(comp$probabilities > 0 & comp$probabilities < 1))
  check("reference_hull_inside", .oe_hull_reference(d$E_c, d$R_c))
  check("reference_hull_outside", !.oe_hull_reference(c(1, 2, 2), c(1, 1, 2)))
  check("reference_hull_collinear", .oe_hull_reference(c(0.5, 1, -1), c(0.5, 1, -1)))

  # Independently difference cumulative logistic probabilities; this routine
  # does not call the extension's stable likelihood or score implementation.
  direct_loglik <- function(p, individual = FALSE) {
    aa <- vapply(z$W, function(w) drop(w %*% p), numeric(n))
    ff <- 1 / (1 + exp(-aa))
    pp <- cbind(ff[, 1], ff[, 2] - ff[, 1], ff[, 3] - ff[, 2], 1 - ff[, 3])
    lli <- log(pp[cbind(seq_len(n), z$y)])
    if (individual) lli else sum(lli)
  }
  check("stable_likelihood_matches_direct", abs(comp$logLik - direct_loglik(theta)) < 1e-9)
  for (h in c(1e-4, 1e-5)) {
    numerical_scores <- numerical_H <- NULL
    numerical_scores <- vapply(seq_along(theta), function(j) {
      dp <- numeric(length(theta)); dp[j] <- h
      (direct_loglik(theta + dp, TRUE) - direct_loglik(theta - dp, TRUE)) / (2 * h)
    }, numeric(n))
    numerical_H <- vapply(seq_along(theta), function(j) {
      dp <- numeric(length(theta)); dp[j] <- h
      -(oe_ppo_components(theta + dp, z)$score - oe_ppo_components(theta - dp, z)$score) / (2 * h)
    }, numeric(length(theta)))
    check(paste0("worker_score_finite_difference_", h), scaled_error(comp$scores, numerical_scores) < 1e-5,
      scaled_error(comp$scores, numerical_scores))
    check(paste0("information_finite_difference_", h), scaled_error(comp$H, numerical_H) < 1e-4,
      scaled_error(comp$H, numerical_H))
  }

  # The package likelihood environment is evaluated at nonoptimal known
  # parameters, making zero-gradient agreement alone insufficient to pass.
  ctl <- ordinal::clm.control(sign.nominal = "positive", sign.location = "negative")
  rho <- ordinal::clm(Y ~ Site + L_c + M_c + B_c + A_c,
    nominal = ~ E_c + R_c, data = d, link = "logit", threshold = "flexible", doFit = FALSE, control = ctl)
  P <- diag(c(rep(1, 3), rep(-1, 6), rep(1, 11)))
  ll_pkg <- rho$clm.nll(rho, par = drop(P %*% theta))
  grad_pkg <- drop(rho$clm.grad(rho)); H_pkg <- rho$clm.hess(rho)
  check("package_nonoptimal_likelihood", abs(ll_pkg + comp$logLik) < 1e-7)
  check("package_nonoptimal_score", scaled_error(comp$score, -drop(P %*% grad_pkg)) < 1e-5)
  check("package_nonoptimal_information", scaled_error(comp$H, P %*% H_pkg %*% P) < 1e-4)

  # Exact PO restriction requires the parameter Jacobian on both derivatives.
  poX <- model.matrix(reformulate(c("Site", .oe_numeric)), d)[, -1L, drop = FALSE]
  po_coef <- setNames(c(rep(0.13, ncol(poX)), -1.1, 0.3, 1.7),
    c(colnames(poX), paste0("threshold:", c("1|2", "2|3", "3|4"))))
  fake_po <- list(valid = TRUE, coef = po_coef)
  restriction <- oe_po_to_ppo(fake_po, z)
  pp <- oe_ppo_components(restriction$theta, z)
  po <- ordinal_components(po_coef, poX, z$y, 4L)
  check("PO_reduction_probability", max(abs(pp$probabilities - po$probabilities)) < 1e-12)
  check("PO_reduction_likelihood", abs(pp$logLik - po$logLik) < 1e-9)
  check("PO_reduction_worker_scores", max(abs(pp$scores %*% restriction$J - po$scores)) < 1e-10)
  check("PO_reduction_information", max(abs(t(restriction$J) %*% pp$H %*% restriction$J - po$H)) < 1e-8)

  # A crossing at an unrealized category is still invalid for that worker.
  dc <- d[1:2, ]; dc$E_c <- c(1, -1); dc$R_c <- 0
  dc$Y <- ordered(c(1, 4), levels = 1:4)
  zc <- oe_ppo_design(dc, FALSE); tc <- setNames(numeric(20), zc$names)
  tc[1:3] <- c(-1, 1, 2); tc[5] <- 4
  check("crossing_unrealized_category_rejected", !is.finite(oe_ppo_components(tc, zc)$logLik))
  tb <- theta; tb[1:9] <- c(-1, -1 + .oe_eps + 5e-7, 1, rep(0, 6))
  boundary <- .oe_assess(tb, z)
  check("boundary_marked_unavailable", !boundary$ok && "boundary_or_near_boundary_unavailable" %in% boundary$reasons)
  huge <- theta; huge["B_c"] <- 10.1
  check("coefficient_growth_flag", "absolute_parameter_above_10_investigation" %in% .oe_assess(huge, z)$reasons)
  ds <- d; ds$B <- 4; ds$B_c <- 1; ds$A <- 3; ds$A_c <- 0
  zs <- oe_ppo_design(ds, FALSE)
  check("staffing_scenarios_leave_ordering_gaps_unchanged",
    max(abs(oe_ppo_components(theta, zs, FALSE)$gaps - comp$gaps)) < 1e-12)

  # Full fits are synthetic and occur only when this test function is called.
  original <- fit_ordinal(d, "Y", .oe_numeric, model_id = "synthetic_original")
  saved$original <- original
  check("synthetic_original_fit_accepted", original$valid, original$status)
  if (isTRUE(original$valid)) {
    source_fixture <- original; source_fixture$model_id <- "q6_main_full"
    source_fixture$analysis_metadata <- list(outcome = "Y", kind = "ordinal", predictors = .oe_numeric,
      row_identifiers = as.character(d$Observation))
    prep <- .prepare_model(d, "Y", .oe_numeric, "Site")
    rejects_source <- function(f) tryCatch({ oe_validate_source_fit(f, prep, "q6_main_full"); FALSE }, error = function(e) TRUE)
    check("matching_saved_source_accepted", !rejects_source(source_fixture))
    wrong <- source_fixture; wrong$model_id <- "q4_main_full"
    check("wrong_saved_model_identifier_rejected", rejects_source(wrong))
    wrong <- source_fixture; wrong$coef <- wrong$coef[c(2, 1, 3:length(wrong$coef))]
    check("swapped_saved_coefficient_order_rejected", rejects_source(wrong))
    wrong <- source_fixture; wrong$analysis_metadata$predictors <- rev(.oe_numeric)
    check("wrong_saved_predictor_metadata_rejected", rejects_source(wrong))
    ppo <- fit_ppo_extension(d, original, "synthetic_ppo")
    saved$ppo <- ppo
    check("PPO_package_constrained_agreement", ppo$valid, ppo$status)
    check("PPO_point_only", identical(ppo$inference_valid, FALSE))
    if (isTRUE(ppo$valid)) {
      check("PPO_recovered_staffing", max(abs(ppo$coef[c("B_c", "A_c")] - theta[c("B_c", "A_c")])) < 0.15)
      pred <- oe_extension_predict(ppo, ds)
      check("PPO_scenario_predictions_coherent", all(pred > 0) && max(abs(rowSums(pred) - 1)) < 1e-8)
    }
    fd <- .oe_factor_data(d)
    fx <- model.matrix(reformulate(c("Site", .oe_factor_predictors)), fd)[, -1L, drop = FALSE]
    check("factor_design_count_and_reference", ncol(fx) == 19L && all(c("B_f1", "B_f2", "B_f4", "B_f5",
      "A_f1", "A_f2", "A_f4", "A_f5") %in% colnames(fx)) && !any(c("B_f3", "A_f3") %in% colnames(fx)))
    reduced <- setNames(numeric(ncol(fx)), colnames(fx))
    shared <- intersect(names(reduced), names(original$coef)); reduced[shared] <- original$coef[shared]
    for (v in c("B", "A")) for (level in c(1, 2, 4, 5))
      reduced[paste0(v, "_f", level)] <- original$coef[paste0(v, "_c")] * (level - 3)
    check("factor_numeric_nesting_design", max(abs(drop(fx %*% reduced) -
      drop(original$design %*% original$coef[colnames(original$design)]))) < 1e-10)
    fac <- fit_ordinal(fd, "Y", .oe_factor_predictors, model_id = "synthetic_factor")
    saved$factor <- fac
    check("factor_fit_accepted", fac$valid, fac$status)
    check("factor_nested_likelihood", isTRUE(fac$valid) && fac$logLik >= original$logLik - 1e-6 * (1 + abs(original$logLik)))
    check("factor_eight_staffing_alternate_agreement", !length(oe_factor_alternate_reasons(fac)))
  }
  factor_names <- paste0(rep(c("B_f", "A_f"), each = 4), rep(c(1, 2, 4, 5), 2))
  regression_fixture <- list(coef = setNames(rep(0, 8), factor_names), alternate = setNames(rep(0, 8), factor_names))
  regression_fixture$alternate["A_f5"] <- 0.005
  check("factor_disagreement_below_generic_tolerance_still_rejected",
    identical(oe_factor_alternate_reasons(regression_fixture), "factor_staffing_alternate_disagreement"))
  regression_fixture$alternate["A_f5"] <- 5e-5
  check("factor_agreement_within_staffing_tolerance_accepted", !length(oe_factor_alternate_reasons(regression_fixture)))
  refused <- tryCatch({ oe_refuse_joint_factor_wald(); FALSE }, error = function(e) grepl("seven-rank", conditionMessage(e)))
  check("eight_restriction_wald_refused", refused)

  # Hand-constructed support counts: 60 workers in all four (3/4,3/4)
  # cells over sites 1:3, seven unsupported workers in site 4, and ten
  # B-only eligible workers in site 5.
  support_rows <- list()
  for (g in 1:3) for (b in 3:4) for (a in 3:4)
    support_rows[[length(support_rows) + 1L]] <- data.frame(Site = g, B = rep(b, 5), A = rep(a, 5))
  support_rows[[length(support_rows) + 1L]] <- data.frame(Site = 4, B = rep(5, 7), A = rep(5, 7))
  for (b in 3:4) support_rows[[length(support_rows) + 1L]] <- data.frame(Site = 5, B = rep(b, 5), A = rep(5, 5))
  sd <- do.call(rbind, support_rows); sd$Observation <- as.character(seq_len(nrow(sd)))
  sd$Y <- ordered(rep(1:4, length.out = nrow(sd)), levels = 1:4)
  for (v in c("L_c", "E_c", "R_c", "M_c")) sd[[v]] <- 0
  sd$B_c <- sd$B - 3; sd$A_c <- sd$A - 3; sd <- .oe_data(sd)
  sup <- oe_staffing_support(sd, 6)
  check("support_hand_count_B", sum(sup$eligible$B) == 70L)
  check("support_hand_count_A", sum(sup$eligible$A) == 60L)
  check("support_intersection_same_population", sum(sup$eligible$AB) == 60L && identical(sup$eligible$AB, sup$eligible$A & sup$eligible$B))
  check("support_site_worker_minima", all(sup$status$supported) && sup$status$sites[sup$status$population == "B"] == 4L)
  check("Q4_has_no_direct_probability_comparison", !"AB" %in% names(oe_staffing_support(sd, 4)$eligible))
  px <- model.matrix(reformulate(c("Site", .oe_numeric)), sd)[, -1L, drop = FALSE]
  b <- setNames(c(rep(0, ncol(px)), -1, 0.3, 1.8), c(colnames(px), paste0("threshold:", c("1|2", "2|3", "3|4"))))
  b[c("B_c", "A_c")] <- c(0.25, 0.4)
  sf <- list(model_id = "supported_test", valid = TRUE, type = "ordinal", status = "accepted", coef = b,
    design = px, predictors = .oe_numeric, n = nrow(sd), G = 5L)
  ps <- oe_probability_summary(list(supported_test = sf), sd, 6, sup)
  check("probability_all_models_accepted", all(ps$status$valid))
  check("probability_four_category_changes_zero", all(vapply(split(ps$contrasts[ps$contrasts$category != "agree_or_strongly_agree", ],
    ps$contrasts$population[ps$contrasts$category != "agree_or_strongly_agree"]), function(x) abs(sum(x$difference)) < 1e-10, logical(1))))
  ab <- ps$contrasts[ps$contrasts$population == "AB", ]
  dd <- sd[sup$eligible$AB, ]
  change <- function(item) {
    d3 <- d4 <- dd; d3[[item]] <- 3; d3[[paste0(item, "_c")]] <- 0
    d4[[item]] <- 4; d4[[paste0(item, "_c")]] <- 1
    colMeans(oe_extension_predict(sf, d4) - oe_extension_predict(sf, d3))
  }
  check("direct_difference_uses_common_population", max(abs(ab$difference[1:4] - (change("A") - change("B")))) < 1e-12)
  check("probability_point_only_columns", !any(grepl("se|lower|upper|p_value|p_reference", names(ps$contrasts))))
  low <- sd[as.character(sd$Site) == "1", ]
  lowps <- oe_probability_summary(list(supported_test = sf), low, 6)
  check("insufficient_support_suppressed", all(!lowps$status$valid) && all(grepl("insufficient", lowps$status$status)))
  failed <- sf; failed$valid <- FALSE; failed$status <- "boundary_or_near_boundary_unavailable"
  bps <- oe_probability_summary(list(failed = failed), sd, 6, sup)
  check("failed_model_predictions_suppressed", all(!bps$status$valid))
  failedpoint <- .oe_point_results(list(model_id = "bad", type = "ppo_ER", valid = FALSE,
    coef = b, n = nrow(sd), G = 5, status = "boundary"), 6)
  check("failed_model_point_estimates_suppressed", all(is.na(failedpoint$estimate)))
  check("Q4_has_no_coefficient_difference", !"A_minus_B" %in% .oe_point_results(sf, 4)$term)
  saved$support <- sup; saved$probabilities <- ps
  ans <- do.call(rbind, records)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(ans, file.path(output_dir, "ordinal_extension_tests.csv"), row.names = FALSE)
    saveRDS(list(checks = ans, objects = saved, seed = seed,
      versions = c(R = as.character(getRversion()), ordinal = as.character(utils::packageVersion("ordinal")))),
      file.path(output_dir, "ordinal_extension_test_objects.rds"))
  }
  if (!all(ans$passed)) stop("Ordinal extension verification failed: ", paste(ans$test[!ans$passed], collapse = "; "))
  ans
}
