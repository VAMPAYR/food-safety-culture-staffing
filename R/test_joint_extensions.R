# Synthetic matrix and saved-object fixtures. No test runs on source.
# Run only after the recorded mathematical and implementation review gates.

.jx_fixture <- function(outcome = "Q6") {
  prefix <- tolower(outcome)
  sites <- rep(as.character(1:8), each = 3L)
  n <- length(sites)
  substantive <- sapply(seq_along(.jx_predictors), function(j) sin(seq_len(n) * j / 7))
  colnames(substantive) <- .jx_predictors
  design <- cbind(sapply(as.character(2:8), function(g) as.numeric(sites == g)), substantive)
  colnames(design) <- c(paste0("Site", 2:8), .jx_predictors)
  param <- c(colnames(design), "threshold:1|2", "threshold:2|3", "threshold:3|4")
  q <- length(param)
  coef <- setNames(rep(0, q), param); coef[c("B_c", "A_c")] <- c(.2, .35)
  coef[grep("^threshold:", names(coef))] <- c(-1, 0, 1)
  H <- diag(q); dimnames(H) <- list(param, param)
  z <- cbind(c(-3, -2, -1, 0, 0, 1, 2, 3), c(-1, 2, -2, 1, -1, 2, -2, 1)) / 10
  S <- matrix(0, 8L, q, dimnames = list(as.character(1:8), param)); S[, .jx_staff] <- z
  scores <- S[sites, , drop = FALSE] / 3; rownames(scores) <- NULL
  V <- (8 / 7) * crossprod(S)
  metadata <- list(outcome = "Y", predictors = .jx_predictors, kind = "ordinal")
  full <- list(model_id = paste0(prefix, "_main_full"), valid = TRUE, type = "ordinal",
    n = n, G = 8L, rows = as.character(seq_len(n)), site = sites,
    outcome = rep(1:4, length.out = n), coef = coef, design = design,
    scores = scores, site_scores = S, H = H, V = V, analysis_metadata = metadata)
  # Nonzero mean deletion displacement distinguishes CV3 from CV3J.
  displacement <- z / 2 + matrix(rep(c(.03, -.02), each = 8L), 8L)
  deleted <- setNames(lapply(as.character(1:8), function(g) {
    keep <- sites != g; kept <- setdiff(as.character(1:8), g)
    site_design <- sapply(kept[-1L], function(s) as.numeric(sites[keep] == s))
    colnames(site_design) <- paste0("Site", kept[-1L])
    X <- cbind(site_design, substantive[keep, , drop = FALSE])
    cf <- c(setNames(rep(0, length(kept) - 1L), colnames(site_design)),
      coef[c(.jx_predictors, grep("^threshold:", names(coef), value = TRUE))])
    cf[.jx_staff] <- cf[.jx_staff] + displacement[as.integer(g), ]
    list(model_id = paste0(prefix, "_omit_site_", g), valid = TRUE, type = "ordinal",
      n = sum(keep), G = 7L, rows = full$rows[keep], site = sites[keep],
      outcome = full$outcome[keep], coef = cf, design = X, analysis_metadata = metadata)
  }), as.character(1:8))
  list(full = full, deletions = deleted, displacement = displacement, site_vectors = z)
}

run_joint_extension_tests <- function(extension_path, output_dir = NULL) {
  source(extension_path, local = environment())
  # Give the fixture the sourced helper constants without modifying globals.
  environment(.jx_fixture) <- environment()
  results <- list()
  check <- function(name, passed, detail = "") {
    results[[length(results) + 1L]] <<- data.frame(test = name, passed = isTRUE(passed),
      detail = as.character(detail), stringsAsFactors = FALSE)
  }
  f <- .jx_fixture()
  ans <- joint_extension_calculations(f$full, f$deletions, "Q6")
  check("fixture_accepted_both_covariances", all(ans$summary$valid), paste(ans$summary$status, collapse = ";"))
  check("all_validation_checks_pass", all(ans$validation$passed))
  check("all_fixed_output_rows_retained", nrow(ans$summary) == 2L && nrow(ans$geometry) == 4L && nrow(ans$contributions) == 16L)
  expected_J <- (7 / 8) * crossprod(f$displacement)
  centered <- sweep(f$displacement, 2, colMeans(f$displacement))
  mean_centered <- (7 / 8) * crossprod(centered)
  correction <- 7 * tcrossprod(colMeans(f$displacement))
  check("full_centering_matches_direct_formula", max(abs(ans$diagnostics$CV3$V - expected_J)) < 1e-12)
  check("full_centering_is_not_mean_centering", max(abs(expected_J - mean_centered)) > 1e-4)
  check("centering_identity", max(abs(expected_J - mean_centered - correction)) < 1e-12)
  expected_S <- (8 / 7) * crossprod(f$site_vectors)
  check("site_score_covariance_direct_formula", max(abs(ans$diagnostics$site$V - expected_S)) < 1e-12)
  b <- unname(f$full$coef[.jx_staff])
  independent_W <- sum(b * solve(expected_J, b))
  j <- ans$summary[ans$summary$covariance_method == "CV3", ]
  check("joint_reference_matches_independent_matrix_calculation", abs(j$W - independent_W) < 1e-12 &&
    abs(j$F_reference - independent_W / 2) < 1e-12 &&
    abs(j$p_reference_uncalibrated - pf(independent_W / 2, 2, 7, lower.tail = FALSE)) < 1e-12)
  check("negative_off_diagonal_contributions_preserved", any(ans$contributions$contribution_BA < 0))
  for (m in c("site", "CV3")) {
    ctab <- ans$contributions[ans$contributions$covariance_method == m, ]
    check(paste0(m, "_shares_sum_one"), abs(sum(ctab$trace_share) - 1) < 1e-12 &&
      abs(sum(ctab$small_eigen_direction_share) - 1) < 1e-12)
    dg <- ans$diagnostics[[m]]
    flipped <- -dg$eigenvectors
    independent_parts <- drop(crossprod(flipped, b))^2 / dg$eigenvalues
    gtab <- ans$geometry[ans$geometry$covariance_method == m, ]
    check(paste0(m, "_eigen_sign_invariance"), max(abs(independent_parts - gtab$W_component)) < 1e-12)
    independent_shares <- vapply(1:8, function(g) drop(crossprod(flipped[, 1L],
      dg$contribution_matrices[, , g] %*% flipped[, 1L])) / dg$eigenvalues[1L], numeric(1))
    check(paste0(m, "_direction_share_sign_invariance"), max(abs(independent_shares - ctab$small_eigen_direction_share)) < 1e-12)
  }
  perm <- f
  order_rows <- rev(seq_along(perm$full$rows))
  for (nm in c("rows", "site", "outcome")) perm$full[[nm]] <- perm$full[[nm]][order_rows]
  perm$full$design <- perm$full$design[order_rows, , drop = FALSE]
  perm$full$scores <- perm$full$scores[order_rows, , drop = FALSE]
  perm$full$site_scores <- perm$full$site_scores[8:1, , drop = FALSE]
  perm$deletions <- perm$deletions[8:1]
  for (g in names(perm$deletions)) {
    ix <- rev(seq_along(perm$deletions[[g]]$rows))
    for (nm in c("rows", "site", "outcome")) perm$deletions[[g]][[nm]] <- perm$deletions[[g]][[nm]][ix]
    perm$deletions[[g]]$design <- perm$deletions[[g]]$design[ix, , drop = FALSE]
  }
  p <- joint_extension_calculations(perm$full, perm$deletions, "Q6")
  check("row_and_site_order_invariance", all(p$summary$valid) && max(abs(p$summary$W - ans$summary$W)) < 1e-12)
  reordered <- f
  ix <- rev(seq_along(reordered$full$coef))
  reordered$full$coef <- reordered$full$coef[ix]
  reordered$full$scores <- reordered$full$scores[, ix, drop = FALSE]
  reordered$full$site_scores <- reordered$full$site_scores[, ix, drop = FALSE]
  reordered$full$H <- reordered$full$H[ix, ix, drop = FALSE]
  reordered$full$V <- reordered$full$V[ix, ix, drop = FALSE]
  for (g in names(reordered$deletions)) reordered$deletions[[g]]$coef <- rev(reordered$deletions[[g]]$coef)
  w <- joint_extension_calculations(reordered$full, reordered$deletions, "Q6")
  check("parameter_order_invariance_with_name_alignment", all(w$summary$valid) &&
    max(abs(w$summary$W - ans$summary$W)) < 1e-12)
  wrong <- f; wrong$deletions[["3"]]$design[1L, "E_c"] <- wrong$deletions[["3"]]$design[1L, "E_c"] + .1
  w <- joint_extension_calculations(wrong$full, wrong$deletions, "Q6")
  check("wrong_deletion_covariates_rejected", !w$summary$valid[2L] && w$summary$valid[1L] &&
    grepl("substantive_design", w$summary$status[2L]))
  wrong <- f; wrong$deletions[["3"]]$design <- cbind(wrong$deletions[["3"]]$design, unexpected = 1)
  w <- joint_extension_calculations(wrong$full, wrong$deletions, "Q6")
  check("extra_deletion_predictor_rejected", !w$summary$valid[2L])
  wrong <- f; wrong$deletions[["3"]]$rows[2L] <- wrong$deletions[["3"]]$rows[1L]
  w <- joint_extension_calculations(wrong$full, wrong$deletions, "Q6")
  check("duplicated_deletion_worker_rejected", !w$summary$valid[2L])
  wrong <- f; wrong$deletions[["2"]]$outcome[1L] <- 4L
  w <- joint_extension_calculations(wrong$full, wrong$deletions, "Q6")
  check("wrong_deletion_outcome_rejected", !w$summary$valid[2L])
  wrong <- f; wrong$deletions[["2"]]$valid <- FALSE
  w <- joint_extension_calculations(wrong$full, wrong$deletions, "Q6")
  check("one_invalid_deletion_invalidates_entire_jackknife", !w$summary$valid[2L] &&
    all(!w$geometry$valid[w$geometry$covariance_method == "CV3"]) &&
    all(!w$contributions$valid[w$contributions$covariance_method == "CV3"]) && w$summary$valid[1L])
  wrong <- f; wrong$deletions <- wrong$deletions[-1L]
  w <- joint_extension_calculations(wrong$full, wrong$deletions, "Q6")
  check("missing_deletion_rejected", !w$summary$valid[2L])
  wrong <- f; names(wrong$deletions)[2L] <- "1"
  w <- joint_extension_calculations(wrong$full, wrong$deletions, "Q6")
  check("duplicate_deletion_rejected", !w$summary$valid[2L])
  wrong <- f; wrong$full$site_scores[1L, "B_c"] <- wrong$full$site_scores[1L, "B_c"] + .1
  w <- joint_extension_calculations(wrong$full, wrong$deletions, "Q6")
  check("saved_site_score_mismatch_rejected", !w$summary$valid[1L] && w$summary$valid[2L])
  wrong <- f; wrong$full$V["B_c", "B_c"] <- NaN
  w <- joint_extension_calculations(wrong$full, wrong$deletions, "Q6")
  check("nonfinite_covariance_rejected", !w$summary$valid[1L])
  wrong <- f
  for (g in as.character(1:8)) wrong$deletions[[g]]$coef[.jx_staff] <-
    wrong$full$coef[.jx_staff] + as.integer(g) / 20 * c(1, 1)
  w <- joint_extension_calculations(wrong$full, wrong$deletions, "Q6")
  check("rank_one_jackknife_has_no_generalized_inverse", !w$summary$valid[2L] && is.na(w$summary$W[2L]))
  equal <- f
  steps <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1), c(1, 0), c(-1, 0), c(0, 1), c(0, -1)) / 10
  for (g in as.character(1:8)) equal$deletions[[g]]$coef[.jx_staff] <-
    equal$full$coef[.jx_staff] + steps[as.integer(g), ]
  w <- joint_extension_calculations(equal$full, equal$deletions, "Q6")
  check("equal_eigenvalues_suppress_directional_shares", w$summary$valid[2L] &&
    all(!w$geometry$orientation_unique[w$geometry$covariance_method == "CV3"]) &&
    all(is.na(w$contributions$small_eigen_direction_share[w$contributions$covariance_method == "CV3"])))
  zero <- f; zero$full$coef[.jx_staff] <- 0
  for (g in as.character(1:8)) zero$deletions[[g]]$coef[.jx_staff] <- f$displacement[as.integer(g), ]
  w <- joint_extension_calculations(zero$full, zero$deletions, "Q6")
  check("zero_coefficient_vector_has_undefined_W_fractions", all(w$summary$valid) &&
    all(w$summary$W == 0) && all(is.na(w$geometry$W_fraction)))
  f4 <- .jx_fixture("Q4")
  a4 <- joint_extension_calculations(f4$full, f4$deletions, "Q4")
  check("Q4_has_only_joint_reference_no_difference_or_Holm", all(a4$summary$valid) && nrow(a4$summary) == 2L &&
    !any(grepl("Holm|A_minus_B|difference", c(names(a4$summary), a4$summary$status))))
  out <- do.call(rbind, results)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(output_dir, "joint_extension_tests.csv")
    if (file.exists(path)) stop("Test output already exists; preserve previous verification.")
    write.csv(out, path, row.names = FALSE)
    saveRDS(list(baseline = ans, Q4 = a4), file.path(output_dir, "joint_extension_test_diagnostics.rds"))
  }
  if (any(!out$passed)) stop("Joint extension synthetic verification failed; inspect returned/saved checks.")
  invisible(out)
}
