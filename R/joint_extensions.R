# Specification 1.4: joint staffing covariance diagnostics.
# Sourcing defines functions only. There are no refits or automatic analyses.
# Run only after the recorded mathematical and implementation review gates.

.jx_staff <- c("B_c", "A_c")
.jx_predictors <- c("L_c", "E_c", "R_c", "M_c", "B_c", "A_c")

.jx_max_error <- function(x, y) {
  if (!identical(dim(x), dim(y)) || length(x) != length(y) ||
      any(!is.finite(x)) || any(!is.finite(y))) return(Inf)
  if (!length(x)) return(0)
  max(abs(x - y))
}

.jx_matrix_tolerance <- function(x) 1e-8 * max(1, abs(x))

.jx_spd <- function(x) {
  failed <- list(valid = FALSE, eigenvalues = c(NA_real_, NA_real_),
                 condition = Inf, symmetric = FALSE)
  if (!is.matrix(x) || !identical(dim(x), c(2L, 2L)) || any(!is.finite(x)))
    return(failed)
  failed$symmetric <- max(abs(x - t(x))) <= .jx_matrix_tolerance(x)
  if (!failed$symmetric) return(failed)
  e <- eigen((x + t(x)) / 2, symmetric = TRUE, only.values = TRUE)$values
  condition <- if (min(e) > 0) max(e) / min(e) else Inf
  list(valid = min(e) > 0 && condition < 1e10,
       eigenvalues = sort(e), condition = condition, symmetric = TRUE)
}

.jx_empty_summary <- function(outcome, method, full, source_run, specification_version) {
  bb <- if (is.numeric(full$coef) && all(.jx_staff %in% names(full$coef)))
    unname(full$coef[.jx_staff]) else c(NA_real_, NA_real_)
  data.frame(outcome = outcome, covariance_method = method,
    full_model_id = if (length(full$model_id) == 1L) full$model_id else NA_character_,
    n = if (length(full$n) == 1L) full$n else NA_integer_,
    G = if (length(full$G) == 1L) full$G else NA_integer_,
    b_B = bb[1L], b_A = bb[2L], var_B = NA_real_, cov_BA = NA_real_,
    var_A = NA_real_, correlation = NA_real_, determinant = NA_real_,
    min_eigenvalue = NA_real_, max_eigenvalue = NA_real_, condition_number = NA_real_,
    W = NA_real_, F_reference = NA_real_, df1 = 2L, df2 = 7L,
    p_reference_uncalibrated = NA_real_, valid = FALSE, status = "unavailable",
    source_run = source_run, specification_version = specification_version,
    stringsAsFactors = FALSE)
}

.jx_empty_geometry <- function(outcome, method) data.frame(
  outcome = outcome, covariance_method = method, eigen_index = 1:2,
  eigenvalue = NA_real_, direction_B = NA_real_, direction_A = NA_real_,
  coefficient_projection = NA_real_, W_component = NA_real_, W_fraction = NA_real_,
  orientation_unique = NA, valid = FALSE, status = "unavailable", stringsAsFactors = FALSE)

.jx_empty_contributions <- function(outcome, method) data.frame(
  outcome = outcome, covariance_method = method, site = as.character(1:8),
  contribution_BB = NA_real_, contribution_BA = NA_real_, contribution_AA = NA_real_,
  trace_share = NA_real_, small_eigen_direction_share = NA_real_,
  displacement_B = NA_real_, displacement_A = NA_real_, orientation_unique = NA,
  valid = FALSE, status = "unavailable", stringsAsFactors = FALSE)

.jx_fit_structure <- function(fit, expected_id, full = NULL) {
  if (!is.list(fit)) return("saved_fit_is_not_a_list")
  reasons <- character()
  if (!isTRUE(fit$valid)) reasons <- c(reasons, "source_fit_not_accepted")
  if (!identical(fit$model_id, expected_id)) reasons <- c(reasons, "unexpected_model_id")
  if (!identical(fit$type, "ordinal")) reasons <- c(reasons, "unexpected_model_type")
  n <- length(fit$rows)
  if (!n || anyNA(fit$rows) || anyDuplicated(fit$rows) ||
      length(fit$site) != n || length(fit$outcome) != n || !identical(as.integer(fit$n), as.integer(n)))
    reasons <- c(reasons, "invalid_rows_sites_or_outcome")
  if (!is.numeric(fit$coef) || !all(.jx_staff %in% names(fit$coef)) ||
      any(!is.finite(fit$coef))) reasons <- c(reasons, "invalid_coefficients")
  kept <- sort(unique(as.character(fit$site)))
  expected_sites <- if (is.null(full)) as.character(1:8) else
    setdiff(as.character(1:8), sub(".*_omit_site_", "", expected_id))
  if (!identical(kept, expected_sites) || !identical(as.integer(fit$G), as.integer(length(expected_sites))))
    reasons <- c(reasons, "unexpected_site_set")
  expected_design <- c(paste0("Site", expected_sites[-1L]), .jx_predictors)
  if (!is.matrix(fit$design) || nrow(fit$design) != n ||
      any(!is.finite(fit$design)) || anyDuplicated(colnames(fit$design)) ||
      !setequal(colnames(fit$design), expected_design)) {
    reasons <- c(reasons, "unexpected_design_columns")
  } else {
    for (label in expected_sites[-1L]) {
      if (!identical(as.numeric(fit$design[, paste0("Site", label)]),
                     as.numeric(as.character(fit$site) == label)))
        reasons <- c(reasons, "site_design_does_not_match_labels")
    }
  }
  # Threshold names in the engine are labels such as threshold:1|2. Check the
  # three distinct threshold names separately from the exact slope design.
  threshold_names <- grep("^threshold:", names(fit$coef), value = TRUE)
  slope_names <- setdiff(names(fit$coef), threshold_names)
  if (length(threshold_names) != 3L || anyDuplicated(names(fit$coef)) ||
      !setequal(slope_names, expected_design)) reasons <- c(reasons, "unexpected_parameter_set")
  meta <- fit$analysis_metadata
  if (!is.list(meta) || !identical(meta$kind, "ordinal") || !identical(meta$outcome, "Y") ||
      !setequal(meta$predictors, .jx_predictors) || length(meta$predictors) != 6L)
    reasons <- c(reasons, "unexpected_analysis_metadata")
  if (!is.null(full)) {
    omitted <- sub(".*_omit_site_", "", expected_id)
    wanted <- full$rows[as.character(full$site) != omitted]
    if (length(fit$rows) != length(wanted) || !setequal(fit$rows, wanted)) {
      reasons <- c(reasons, "deletion_rows_do_not_match")
    } else {
      idx <- match(fit$rows, full$rows)
      if (!identical(fit$outcome, full$outcome[idx])) reasons <- c(reasons, "deletion_outcome_does_not_match")
      if (!identical(as.character(fit$site), as.character(full$site[idx])))
        reasons <- c(reasons, "deletion_site_labels_do_not_match")
      if (is.matrix(fit$design) && is.matrix(full$design) &&
          all(.jx_predictors %in% colnames(fit$design)) && all(.jx_predictors %in% colnames(full$design))) {
        if (!isTRUE(all.equal(unname(fit$design[, .jx_predictors, drop = FALSE]),
                             unname(full$design[idx, .jx_predictors, drop = FALSE]), tolerance = 0)))
          reasons <- c(reasons, "deletion_substantive_design_does_not_match")
      } else reasons <- c(reasons, "deletion_substantive_design_unavailable")
    }
  }
  unique(reasons)
}

joint_extension_calculations <- function(full, deletions, outcome,
                                         source_run = "in_memory", specification_version = "1.4") {
  if (!outcome %in% c("Q6", "Q4")) stop("outcome must be Q6 or Q4")
  if (!is.list(full)) full <- list(valid = FALSE, status = "saved_full_fit_is_not_a_list")
  prefix <- tolower(outcome)
  checks <- list()
  add_check <- function(method, name, value = NA_real_, tolerance = NA_real_, passed, reason = "") {
    checks[[length(checks) + 1L]] <<- data.frame(outcome = outcome,
      method = method, check_name = name, observed_error_or_value = value,
      tolerance = tolerance, passed = isTRUE(passed), reason = reason, stringsAsFactors = FALSE)
    isTRUE(passed)
  }
  summaries <- geometries <- contributions <- diagnostics <- list()
  full_reasons <- .jx_fit_structure(full, paste0(prefix, "_main_full"))
  add_check("both", "full_fit_structure", passed = !length(full_reasons), reason = paste(full_reasons, collapse = ";"))
  bb <- if (!length(full_reasons)) unname(full$coef[.jx_staff]) else rep(NA_real_, 2L)

  for (method in c("site", "CV3")) {
    summary <- .jx_empty_summary(outcome, method, full, source_run, specification_version)
    geometry <- .jx_empty_geometry(outcome, method)
    contribution <- .jx_empty_contributions(outcome, method)
    failure <- full_reasons
    V <- matrix(NA_real_, 2L, 2L, dimnames = list(.jx_staff, .jx_staff))
    pieces <- array(NA_real_, c(2L, 2L, 8L), dimnames = list(.jx_staff, .jx_staff, as.character(1:8)))
    d <- matrix(NA_real_, 8L, 2L, dimnames = list(as.character(1:8), .jx_staff))
    site_score_error <- NA_real_
    if (!length(failure) && method == "site") {
      param <- names(full$coef)
      q <- length(param)
      arrays_ok <- is.matrix(full$scores) && identical(dim(full$scores), c(length(full$rows), q)) &&
        identical(colnames(full$scores), param) && all(is.finite(full$scores)) &&
        is.matrix(full$H) && identical(dim(full$H), c(q, q)) &&
        identical(rownames(full$H), param) && identical(colnames(full$H), param) && all(is.finite(full$H)) &&
        is.matrix(full$site_scores) && identical(dim(full$site_scores), c(8L, q)) &&
        setequal(rownames(full$site_scores), as.character(1:8)) &&
        identical(colnames(full$site_scores), param) && all(is.finite(full$site_scores)) &&
        is.matrix(full$V) && identical(dim(full$V), c(q, q)) &&
        identical(rownames(full$V), param) && identical(colnames(full$V), param) && all(is.finite(full$V))
      add_check(method, "score_information_array_names_and_dimensions", passed = arrays_ok)
      if (!arrays_ok) failure <- c(failure, "score_information_arrays_invalid") else {
        S <- rowsum(full$scores, as.character(full$site), reorder = TRUE)
        S <- S[as.character(1:8), , drop = FALSE]
        site_score_error <- .jx_max_error(S, full$site_scores[as.character(1:8), , drop = FALSE])
        score_ok <- add_check(method, "site_score_aggregation", site_score_error,
          .jx_matrix_tolerance(S), site_score_error <= .jx_matrix_tolerance(S))
        stationarity <- max(abs(colSums(full$scores))) / full$n
        stationary <- add_check(method, "accepted_score_tolerance", stationarity, 1e-6, stationarity < 1e-6)
        H_symmetric <- max(abs(full$H - t(full$H))) <= .jx_matrix_tolerance(full$H)
        H_chol <- if (H_symmetric) tryCatch(chol(full$H), error = function(e) NULL) else NULL
        H_eig <- if (!is.null(H_chol)) eigen(full$H, symmetric = TRUE, only.values = TRUE)$values else c(0, 1)
        H_condition <- if (min(H_eig) > 0) max(H_eig) / min(H_eig) else Inf
        H_ok <- add_check(method, "information_positive_definite_and_condition", H_condition, 1e10,
          !is.null(H_chol) && H_condition <= 1e10)
        if (!score_ok) failure <- c(failure, "site_score_aggregation_mismatch")
        if (!stationary) failure <- c(failure, "accepted_score_tolerance_failed")
        if (!H_ok) failure <- c(failure, "information_invalid")
        if (!length(failure)) {
          # Each column solves H x = S_g. Include every nuisance parameter in
          # the solve, then select the two staffing entries.
          solved <- backsolve(H_chol, forwardsolve(t(H_chol), t(S)))
          u <- solved[match(.jx_staff, param), , drop = FALSE]
          for (g in 1:8) pieces[, , g] <- (8 / 7) * tcrossprod(u[, g])
          V <- full$V[.jx_staff, .jx_staff, drop = FALSE]
        }
      }
    }
    if (!length(failure) && method == "CV3") {
      number_ok <- length(deletions) == 8L && !is.null(names(deletions)) &&
        !anyDuplicated(names(deletions)) && setequal(names(deletions), as.character(1:8))
      add_check(method, "eight_distinct_named_deletions", length(deletions), 8, number_ok)
      if (!number_ok) failure <- c(failure, "exactly_eight_named_deletions_required") else {
        deletion_reasons <- character()
        for (g in as.character(1:8)) {
          f <- deletions[[g]]
          reasons <- .jx_fit_structure(f, paste0(prefix, "_omit_site_", g), full)
          add_check(method, paste0("deletion_", g, "_fit_rows_outcome_design"),
                    passed = !length(reasons), reason = paste(reasons, collapse = ";"))
          if (length(reasons)) deletion_reasons <- c(deletion_reasons, paste0("site_", g, ":", reasons))
          else d[g, ] <- unname(f$coef[.jx_staff]) - bb
        }
        add_check(method, "valid_deletion_count", sum(complete.cases(d)), 8, all(complete.cases(d)))
        failure <- c(failure, deletion_reasons)
        if (!length(failure)) {
          for (g in 1:8) pieces[, , g] <- (7 / 8) * tcrossprod(d[g, ])
          V <- (7 / 8) * crossprod(d)
          dimnames(V) <- list(.jx_staff, .jx_staff)
          contribution$displacement_B <- d[, 1L]
          contribution$displacement_A <- d[, 2L]
        }
      }
    }

    if (!length(failure)) {
      matrix_status <- .jx_spd(V)
      cov_ok <- add_check(method, "covariance_positive_definite_and_condition", matrix_status$condition,
        1e10, matrix_status$valid)
      add_check(method, "covariance_symmetry", if (all(is.finite(V))) max(abs(V - t(V))) else Inf,
        if (all(is.finite(V))) .jx_matrix_tolerance(V) else NA_real_, matrix_status$symmetric)
      reconstruction <- apply(pieces, c(1L, 2L), sum)
      reconstruction_error <- .jx_max_error(V, reconstruction)
      recon_ok <- add_check(method, "covariance_contribution_reconstruction", reconstruction_error,
        .jx_matrix_tolerance(V), reconstruction_error <= .jx_matrix_tolerance(V))
      if (!cov_ok) failure <- c(failure, "covariance_not_positive_definite_or_ill_conditioned")
      if (!recon_ok) failure <- c(failure, "covariance_reconstruction_failed")
    }
    if (!length(failure)) {
      eig <- eigen((V + t(V)) / 2, symmetric = TRUE)
      order_eig <- order(eig$values)
      lambda <- eig$values[order_eig]
      vectors <- eig$vectors[, order_eig, drop = FALSE]
      for (r in 1:2) {
        first <- which(abs(vectors[, r]) > 1e-12)[1L]
        if (vectors[first, r] < 0) vectors[, r] <- -vectors[, r]
      }
      orientation_unique <- (lambda[2L] - lambda[1L]) > 1e-8 * lambda[2L]
      projection <- drop(crossprod(vectors, bb))
      parts <- projection^2 / lambda
      R <- chol(V)
      solved_b <- backsolve(R, forwardsolve(t(R), bb))
      W <- sum(bb * solved_b)
      W_error <- abs(W - sum(parts))
      orth_error <- max(abs(crossprod(vectors) - diag(2L)))
      eig_error <- .jx_max_error(V, vectors %*% diag(lambda, 2L) %*% t(vectors))
      raw_trace <- vapply(1:8, function(g) sum(diag(pieces[, , g])) / sum(diag(V)), numeric(1))
      raw_direction <- if (orientation_unique) vapply(1:8, function(g)
        drop(crossprod(vectors[, 1L], pieces[, , g] %*% vectors[, 1L])) / lambda[1L], numeric(1)) else rep(NA_real_, 8L)
      numeric_checks <- c(
        add_check(method, "eigen_reconstruction", eig_error, .jx_matrix_tolerance(V), eig_error <= .jx_matrix_tolerance(V)),
        add_check(method, "eigenvector_orthogonality", orth_error, 1e-8, orth_error <= 1e-8),
        add_check(method, "W_reconstruction", W_error, 1e-8 * max(1, abs(W)), W_error <= 1e-8 * max(1, abs(W))),
        add_check(method, "trace_share_sum", abs(sum(raw_trace) - 1), 1e-8, abs(sum(raw_trace) - 1) <= 1e-8),
        add_check(method, "trace_share_nonnegative", min(raw_trace), -1e-12, min(raw_trace) >= -1e-12))
      if (orientation_unique) numeric_checks <- c(numeric_checks,
        add_check(method, "small_direction_share_sum", abs(sum(raw_direction) - 1), 1e-8,
          abs(sum(raw_direction) - 1) <= 1e-8),
        add_check(method, "small_direction_share_nonnegative", min(raw_direction), -1e-12,
          min(raw_direction) >= -1e-12))
      else add_check(method, "nonunique_orientation_suppressed", lambda[2L] - lambda[1L],
        1e-8 * lambda[2L], all(is.na(raw_direction)), "Method-specific orientation is nonunique; no directional shares or interpretation.")
      if (!all(numeric_checks)) failure <- c(failure, "derived_covariance_validation_failed") else {
        summary$var_B <- V[1L, 1L]; summary$cov_BA <- V[1L, 2L]; summary$var_A <- V[2L, 2L]
        summary$correlation <- V[1L, 2L] / sqrt(V[1L, 1L] * V[2L, 2L])
        summary$determinant <- prod(lambda)
        summary$min_eigenvalue <- lambda[1L]; summary$max_eigenvalue <- lambda[2L]
        summary$condition_number <- lambda[2L] / lambda[1L]
        summary$W <- W; summary$F_reference <- W / 2
        summary$p_reference_uncalibrated <- pf(W / 2, 2, 7, lower.tail = FALSE)
        summary$valid <- TRUE; summary$status <- "accepted_uncalibrated_F2_7_reference_sensitivity"
        geometry$eigenvalue <- lambda; geometry$direction_B <- vectors[1L, ]; geometry$direction_A <- vectors[2L, ]
        geometry$coefficient_projection <- projection; geometry$W_component <- parts
        geometry$W_fraction <- if (W > 1e-12) parts / W else rep(NA_real_, 2L)
        geometry$orientation_unique <- orientation_unique; geometry$valid <- TRUE
        geometry$status <- if (orientation_unique) "accepted_diagnostic_directions" else "accepted_orientation_nonunique_no_directional_interpretation"
        contribution$contribution_BB <- pieces[1L, 1L, ]; contribution$contribution_BA <- pieces[1L, 2L, ]
        contribution$contribution_AA <- pieces[2L, 2L, ]; contribution$trace_share <- pmax(0, raw_trace)
        contribution$small_eigen_direction_share <- if (orientation_unique) pmax(0, raw_direction) else NA_real_
        contribution$orientation_unique <- orientation_unique; contribution$valid <- TRUE
        contribution$status <- if (orientation_unique) "accepted_method_specific_covariance_shares" else "accepted_trace_shares_direction_nonunique"
        diagnostics[[method]] <- list(V = V, contribution_matrices = pieces,
          displacements = if (method == "CV3") d else NULL, eigenvalues = lambda,
          eigenvectors = vectors, raw_trace_share = raw_trace,
          raw_small_direction_share = raw_direction, site_score_error = site_score_error)
      }
    }
    if (length(failure)) {
      status <- paste(unique(failure), collapse = ";")
      summary$status <- geometry$status <- contribution$status <- status
      diagnostics[[method]] <- list(V = V, contribution_matrices = pieces,
        displacements = d, failure = unique(failure), site_score_error = site_score_error)
    }
    summaries[[method]] <- summary; geometries[[method]] <- geometry; contributions[[method]] <- contribution
  }
  list(summary = do.call(rbind, summaries), geometry = do.call(rbind, geometries),
       contributions = do.call(rbind, contributions), validation = do.call(rbind, checks),
       diagnostics = diagnostics)
}

run_joint_extension <- function(source_run, outdir, package_dir = dirname(dirname(source_run)),
                                specification_version = "1.4") {
  source_run <- normalizePath(source_run, winslash = "/", mustWork = TRUE)
  package_dir <- normalizePath(package_dir, winslash = "/", mustWork = TRUE)
  if (!dir.exists(outdir) && !dir.create(outdir, recursive = TRUE)) stop("Cannot create output directory.")
  outdir <- normalizePath(outdir, winslash = "/", mustWork = TRUE)
  target_names <- c("joint_covariance_summary.csv", "joint_covariance_geometry.csv",
    "jurisdiction_covariance_contributions.csv", "joint_covariance_validation.csv",
    "joint_extension_diagnostics.rds", "joint_extension_source_files.csv", "joint_extension_provenance.rds")
  provenance_dir <- file.path(outdir, "joint_provenance")
  if (any(file.exists(file.path(outdir, target_names))) || dir.exists(provenance_dir))
    stop("Joint extension outputs already exist; do not overwrite a historical run.")
  code_dir <- file.path(package_dir, "R")
  source_files <- file.path(code_dir, c("model_engine.R", "run_analysis.R", "joint_extensions.R", "test_joint_extensions.R"))
  spec_file <- file.path(package_dir, "specification", "mathematical_specification", "main.tex")
  if (!all(file.exists(source_files))) stop("Required implementation files are missing.")
  helpers <- new.env(parent = environment())
  source(source_files[1L], local = helpers)
  source(source_files[2L], local = helpers)
  fit_ids <- unlist(lapply(c("q6", "q4"), function(prefix)
    c(paste0(prefix, "_main_full"), paste0(prefix, "_omit_site_", 1:8))), use.names = FALSE)
  fit_paths <- file.path(source_run, "fits", paste0(fit_ids, ".rds"))
  provenance <- data.frame(model_id = fit_ids, file = fit_paths, exists = file.exists(fit_paths),
    readable = FALSE, read_status = "not_attempted", stringsAsFactors = FALSE)
  fits <- setNames(vector("list", length(fit_ids)), fit_ids)
  for (j in seq_along(fit_paths)) {
    f <- tryCatch(readRDS(fit_paths[j]), error = function(e) e)
    if (inherits(f, "error")) {
      provenance$read_status[j] <- conditionMessage(f)
      fits[[j]] <- list(model_id = fit_ids[j], valid = FALSE, status = conditionMessage(f))
    } else if (!is.list(f)) {
      provenance$read_status[j] <- "RDS object is not a saved fit list"
      fits[[j]] <- list(model_id = fit_ids[j], valid = FALSE, status = provenance$read_status[j])
    } else {
      provenance$readable[j] <- TRUE; provenance$read_status[j] <- "read"
      fits[[j]] <- f
    }
  }
  results <- setNames(lapply(c("Q6", "Q4"), function(outcome) {
    prefix <- tolower(outcome)
    ds <- fits[paste0(prefix, "_omit_site_", 1:8)]; names(ds) <- as.character(1:8)
    joint_extension_calculations(fits[[paste0(prefix, "_main_full")]], ds, outcome,
      source_run = source_run, specification_version = specification_version)
  }), c("Q6", "Q4"))
  table_names <- c(summary = "joint_covariance_summary.csv", geometry = "joint_covariance_geometry.csv",
    contributions = "jurisdiction_covariance_contributions.csv", validation = "joint_covariance_validation.csv")
  for (name in names(table_names)) write.csv(do.call(rbind, lapply(results, `[[`, name)),
    file.path(outdir, table_names[[name]]), row.names = FALSE, na = "")
  saveRDS(lapply(results, `[[`, "diagnostics"), file.path(outdir, "joint_extension_diagnostics.rds"))
  write.csv(provenance, file.path(outdir, "joint_extension_source_files.csv"), row.names = FALSE)
  dir.create(provenance_dir)
  # The internal specification is intentionally absent from the public code
  # release. Bind its hash when available; preserve the version either way.
  extra_files <- c(spec_file, file.path(source_run, c("provenance_sha256.csv", "provenance_md5.csv")))
  helpers$write_sha256_manifest(c(fit_paths[file.exists(fit_paths)], source_files,
    extra_files[file.exists(extra_files)]), provenance_dir)
  saveRDS(list(source_run = source_run, package_dir = package_dir, specification_version = specification_version,
    timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), source_files = provenance,
    internal_specification_available = file.exists(spec_file),
    session_info = sessionInfo(), statement = "Saved fits reused; no refits. All F probabilities are uncalibrated references."),
    file.path(outdir, "joint_extension_provenance.rds"))
  invisible(results)
}
