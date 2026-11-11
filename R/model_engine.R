# Specification 1.1/1.2: cumulative-logit model and site-based inference.
# Dependencies: R base/stats and the recommended MASS package.
# No analysis executes when this file is sourced.

.capture <- function(expr) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(value = if (inherits(value, "error")) NULL else value,
       error = if (inherits(value, "error")) conditionMessage(value) else NULL,
       warnings = unique(warnings))
}

.log_expit <- function(x) -pmax(-x, 0) - log1p(exp(-abs(x)))
.log_one_minus_exp <- function(x) {
  # x must be strictly negative; avoid cancellation near zero.
  ans <- numeric(length(x)); small <- x < -log(2)
  ans[small] <- log1p(-exp(x[small]))
  ans[!small] <- log(-expm1(x[!small])); ans
}

ordinal_components <- function(theta, X, y, K, derivatives = TRUE) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  stopifnot(length(theta) == p + K - 1L, length(y) == n,
            all(y %in% seq_len(K)))
  tau <- theta[p + seq_len(K - 1L)]
  if (any(!is.finite(theta)) || any(diff(tau) <= 0)) {
    return(list(logLik = -Inf, probabilities = matrix(NA_real_, n, K)))
  }
  eta <- drop(X %*% theta[seq_len(p)])
  logp <- matrix(NA_real_, n, K)
  logp[, 1L] <- .log_expit(tau[1L] - eta)
  logp[, K] <- .log_expit(eta - tau[K - 1L])
  if (K > 2L) for (k in 2L:(K - 1L)) {
    a <- tau[k - 1L] - eta; b <- tau[k] - eta
    logp[, k] <- .log_expit(b) + .log_expit(-a) +
      .log_one_minus_exp(rep(tau[k - 1L] - tau[k], n))
  }
  ll_i <- logp[cbind(seq_len(n), y)]
  result <- list(logLik = sum(ll_i), logLik_i = ll_i,
                 probabilities = exp(logp), eta = eta)
  if (!derivatives) return(result)
  q <- length(theta); scores <- matrix(0, n, q)
  H <- matrix(0, q, q)
  # Natural parameters are slopes followed by ordered thresholds. These
  # derivatives use stable logistic identities rather than division by p(y).
  for (i in seq_len(n)) {
    k <- y[i]; lo <- k - 1L; hi <- k
    Fa <- if (lo > 0L) plogis(tau[lo] - eta[i]) else 0
    Fb <- if (hi < K) plogis(tau[hi] - eta[i]) else 1
    fa <- Fa * (1 - Fa); fb <- Fb * (1 - Fb)
    r <- if (lo > 0L && hi < K) 1 / expm1(tau[hi] - tau[lo]) else 0
    rr <- r * (1 + r)
    scores[i, seq_len(p)] <- X[i, ] * (Fa + Fb - 1)
    H[seq_len(p), seq_len(p)] <- H[seq_len(p), seq_len(p), drop = FALSE] +
      tcrossprod(X[i, ]) * (fa + fb)
    if (lo > 0L) {
      j <- p + lo
      scores[i, j] <- -Fa - r
      H[j, j] <- H[j, j] + fa + rr
      H[seq_len(p), j] <- H[seq_len(p), j] - X[i, ] * fa
      H[j, seq_len(p)] <- H[seq_len(p), j]
    }
    if (hi < K) {
      j <- p + hi
      scores[i, j] <- 1 - Fb + r
      H[j, j] <- H[j, j] + fb + rr
      H[seq_len(p), j] <- H[seq_len(p), j] - X[i, ] * fb
      H[j, seq_len(p)] <- H[seq_len(p), j]
    }
    if (lo > 0L && hi < K) {
      H[p + lo, p + hi] <- H[p + lo, p + hi] - rr
      H[p + hi, p + lo] <- H[p + lo, p + hi]
    }
  }
  if (!is.null(names(theta))) {
    colnames(scores) <- names(theta); dimnames(H) <- list(names(theta), names(theta))
  }
  result$scores <- scores; result$score <- colSums(scores); result$H <- H
  result
}

.to_natural <- function(phi, p) {
  t <- phi[-seq_len(p)]
  c(phi[seq_len(p)], cumsum(c(t[1L], exp(t[-1L]))))
}
.to_unconstrained <- function(theta, p) {
  t <- theta[-seq_len(p)]
  c(theta[seq_len(p)], t[1L], log(diff(t)))
}
.phi_gradient <- function(score, theta, p) {
  st <- score[-seq_len(p)]; tau <- theta[-seq_len(p)]
  tail_sums <- rev(cumsum(rev(st)))
  c(score[seq_len(p)], sum(st), diff(tau) * tail_sums[-1L])
}

.matrix_check <- function(A, max_condition = Inf) {
  if (!is.matrix(A) || nrow(A) != ncol(A) || any(!is.finite(A)))
    return(list(valid = FALSE, rank = NA_integer_, condition = Inf, min_eigen = NA_real_))
  e <- eigen((A + t(A)) / 2, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(abs(e)); rank <- if (scale == 0) 0L else sum(abs(e) > scale * 1e-10)
  condition <- if (min(e) > 0) max(e) / min(e) else Inf
  list(valid = min(e) > 0 && condition < max_condition,
       rank = rank, condition = condition, min_eigen = min(e))
}

site_sandwich <- function(scores, H, site) {
  # H and scores are sums on the worker likelihood scale. The only finite
  # cluster multiplier is G/(G-1); there is no n/(n-p) adjustment.
  scores <- as.matrix(scores)
  S <- rowsum(scores, as.character(site), reorder = TRUE)
  G <- nrow(S); info <- .matrix_check(H)
  if (G < 2L || !info$valid) return(list(valid = FALSE, V = NULL,
      site_scores = S, G = G, information = info))
  Hinv <- chol2inv(chol(H))
  V <- G / (G - 1) * Hinv %*% crossprod(S) %*% Hinv
  V <- (V + t(V)) / 2; dimnames(V) <- dimnames(H)
  list(valid = all(is.finite(V)), V = V, site_scores = S, G = G,
       information = info, rank = .matrix_check(V)$rank,
       score_sum = colSums(S))
}

.prepare_model <- function(data, outcome, predictors, site, binary = FALSE) {
  required <- unique(c(outcome, predictors, site))
  if (!all(required %in% names(data))) stop("Missing model columns: ",
      paste(setdiff(required, names(data)), collapse = ", "))
  if (anyNA(data[, required, drop = FALSE])) stop("Model input contains missing values; define rows explicitly.")
  d <- data; sites <- unique(as.character(d[[site]]))
  numeric_sites <- suppressWarnings(as.numeric(sites))
  sites <- if (all(is.finite(numeric_sites))) sites[order(numeric_sites)] else sort(sites)
  d[[site]] <- factor(as.character(d[[site]]), levels = sites)
  if (length(sites) < 2L) stop("At least two observed sites are required.")
  predictors <- setdiff(unique(predictors), site)
  # Treatment contrasts are explicit and unaffected by session options.
  for (nm in c(site, predictors)) if (is.factor(d[[nm]])) {
    d[[nm]] <- droplevels(d[[nm]])
    if (nlevels(d[[nm]]) < 2L) stop("Constant factor: ", nm)
    contrasts(d[[nm]]) <- contr.treatment(levels(d[[nm]]), base = 1L)
  }
  formula <- reformulate(c(site, predictors), response = outcome)
  if (binary) {
    yy <- d[[outcome]]
    if (is.factor(yy)) yy <- as.character(yy)
    yy <- suppressWarnings(as.numeric(yy))
    if (!all(yy %in% c(0, 1)) || length(unique(yy)) != 2L)
      stop("Binary outcome must contain both 0 and 1.")
    d[[outcome]] <- yy; K <- 2L
  } else {
    yy <- d[[outcome]]
    lev <- if (is.factor(yy)) levels(yy) else sort(unique(yy))
    if (length(lev) < 3L) stop("Ordinal outcome needs at least three prespecified levels.")
    d[[outcome]] <- ordered(yy, levels = lev)
    if (any(table(d[[outcome]]) == 0L)) stop("An ordinal outcome level is empty.")
    K <- length(lev); yy <- as.integer(d[[outcome]])
  }
  Xfull <- model.matrix(formula, d)
  if (any(!is.finite(Xfull))) stop("Nonfinite design matrix.")
  if (qr(Xfull)$rank != ncol(Xfull)) stop("Design matrix is not full column rank.")
  X <- if (binary) Xfull else Xfull[, colnames(Xfull) != "(Intercept)", drop = FALSE]
  list(data = d, formula = formula, X = X, y = yy, K = K,
       site = d[[site]], sites = sites,
       rows = if ("Observation" %in% names(d)) as.character(d$Observation) else rownames(d))
}

.failed_fit <- function(model_id, status, n = NA_integer_, diagnostics = list()) {
  structure(list(model_id = model_id, coef = numeric(), coefficients = numeric(),
    theta = numeric(), V = NULL, vcov_site = NULL, valid = FALSE,
    inference_valid = FALSE, status = status, logLik = NA_real_, n = n,
    G = NA_integer_, probabilities = NULL, parameter_names = character(),
    diagnostics = diagnostics, rows = character()), class = "site_model")
}

.finish_fit <- function(prep, theta, comp, model_id, attempts, reasons,
                        investigation, type, alternate) {
  sandwich <- site_sandwich(comp$scores, comp$H, prep$site)
  if (!sandwich$valid) reasons <- c(reasons, "site_covariance_unavailable")
  if (length(investigation)) reasons <- c(reasons, "investigation_required")
  valid <- !length(reasons)
  all_sites <- identical(sort(as.character(prep$sites)), as.character(1:8))
  status <- if (valid) "accepted" else paste(unique(reasons), collapse = ";")
  structure(list(model_id = model_id, coef = theta, theta = theta,
    coefficients = theta, V = sandwich$V, vcov_site = sandwich$V,
    valid = valid, inference_valid = valid && all_sites && sandwich$valid,
    status = status, logLik = comp$logLik, n = nrow(prep$X), G = length(prep$sites),
    probabilities = comp$probabilities, parameter_names = names(theta),
    scores = comp$scores, H = comp$H, site_scores = sandwich$site_scores,
    design = prep$X, outcome = prep$y, rows = prep$rows, site = prep$site,
    type = type, alternate = alternate,
    diagnostics = list(attempts = attempts, reasons = unique(reasons),
      investigation = investigation, max_score_per_n = max(abs(comp$score)) / nrow(prep$X),
      information = sandwich$information, sandwich_rank = sandwich$rank,
      expected_rank_upper_bound = length(prep$sites) - 1L,
      all_eight_sites = all_sites,
      inference_status = if (!all_sites) "All eight sites are required for reference calculations." else status)),
    class = "site_model")
}

.refine_ordinal_primary <- function(z, theta, convergence, attempts) {
  p <- ncol(z$X); n <- nrow(z$X)
  comp <- ordinal_components(theta, z$X, z$y, z$K)
  score <- if (is.null(comp$score)) Inf else max(abs(comp$score)) / n
  attempts$primary$logLik <- comp$logLik
  attempts$primary$max_score_per_n <- score
  attempts$primary$theta <- theta
  attempts$selected_primary <- "primary"
  info <- if (is.null(comp$H)) list(valid = FALSE, condition = Inf) else .matrix_check(comp$H)
  probability_ok <- all(is.finite(comp$probabilities)) &&
    all(comp$probabilities >= 0 & comp$probabilities <= 1) &&
    max(abs(rowSums(comp$probabilities) - 1)) <= 1e-8
  # A generic numerical restart for a converged finite candidate whose only
  # primary numerical defect is the score criterion. No data or model changes.
  eligible <- isTRUE(convergence == 0) && all(is.finite(theta)) &&
    is.finite(comp$logLik) && is.finite(score) && score >= 1e-6 &&
    info$valid && info$condition <= 1e10 && max(abs(theta)) <= 10 &&
    min(diff(theta[p + seq_len(z$K - 1L)])) >= 1e-6 && probability_ok
  if (eligible) {
    retry <- .capture(MASS::polr(z$formula, data = z$data, method = "logistic",
      start = unname(theta), Hess = FALSE, model = TRUE,
      control = list(maxit = 10000L, reltol = 1e-13)))
    record <- list(solver = "MASS::polr BFGS score refinement", trigger = "primary_score_tolerance",
      start = theta, control = list(maxit = 10000L, reltol = 1e-13),
      warnings = retry$warnings, error = retry$error, selected = FALSE)
    record$allowed_logLik_decrease <- 100 * .Machine$double.eps * (1 + abs(comp$logLik))
    if (!is.null(retry$value)) {
      candidate <- retry$value
      candidate_theta <- c(coef(candidate)[colnames(z$X)],
        setNames(candidate$zeta, paste0("threshold:", names(candidate$zeta))))
      cc <- ordinal_components(candidate_theta, z$X, z$y, z$K)
      candidate_score <- if (is.null(cc$score)) Inf else max(abs(cc$score)) / n
      record$theta <- candidate_theta; record$convergence <- candidate$convergence
      record$logLik <- cc$logLik; record$max_score_per_n <- candidate_score
      record$logLik_change <- cc$logLik - comp$logLik
      if (isTRUE(candidate$convergence == 0) && all(is.finite(candidate_theta)) &&
          is.finite(cc$logLik) && is.finite(candidate_score) && candidate_score < score &&
          cc$logLik >= comp$logLik - record$allowed_logLik_decrease) {
        theta <- candidate_theta; comp <- cc; convergence <- candidate$convergence
        record$selected <- TRUE; attempts$selected_primary <- "primary_refinement"
      }
    }
    attempts$primary_refinement <- record
  }
  list(theta = theta, comp = comp, convergence = convergence, attempts = attempts)
}

fit_ordinal <- function(data, outcome, predictors, site = "Site", model_id = "ordinal") {
  prepared <- .capture(.prepare_model(data, outcome, predictors, site))
  if (is.null(prepared$value)) return(.failed_fit(model_id, prepared$error, nrow(data)))
  z <- prepared$value; p <- ncol(z$X); n <- nrow(z$X); K <- z$K
  primary <- .capture(MASS::polr(z$formula, data = z$data, method = "logistic",
    Hess = FALSE, model = TRUE, control = list(maxit = 5000L, reltol = 1e-10)))
  attempts <- list(primary = list(solver = "MASS::polr BFGS", warnings = primary$warnings,
    error = primary$error))
  if (is.null(primary$value)) return(.failed_fit(model_id, "primary_fit_failed", n,
      list(attempts = attempts)))
  m <- primary$value
  theta <- c(coef(m)[colnames(z$X)], setNames(m$zeta, paste0("threshold:", names(m$zeta))))
  attempts$primary$convergence <- m$convergence
  attempts$primary$theta <- theta
  comp <- ordinal_components(theta, z$X, z$y, K)
  if (!is.finite(comp$logLik) || is.null(comp$score))
    return(.failed_fit(model_id, "nonfinite_primary_parameters_or_likelihood", n, list(attempts = attempts)))
  refinement <- .refine_ordinal_primary(z, theta, m$convergence, attempts)
  theta <- refinement$theta; comp <- refinement$comp; attempts <- refinement$attempts
  primary_convergence <- refinement$convergence
  # Independent implementation and perturbed start. Tighter tolerance also
  # tests parameter growth when apparent convergence masks separation.
  phi0 <- .to_unconstrained(theta, p) + 0.03 * sin(seq_along(theta))
  alternate <- .capture(optim(phi0, function(phi) {
    -ordinal_components(.to_natural(phi, p), z$X, z$y, K, FALSE)$logLik
  }, gr = function(phi) {
    natural <- .to_natural(phi, p)
    cc <- ordinal_components(natural, z$X, z$y, K)
    -.phi_gradient(cc$score, natural, p)
  }, method = "BFGS", control = list(maxit = 10000L, reltol = 1e-12)))
  attempts$alternate <- list(solver = "independent likelihood BFGS", warnings = alternate$warnings,
    error = alternate$error)
  reasons <- character(); investigation <- character(); alt_theta <- NULL
  if (!isTRUE(primary_convergence == 0)) reasons <- c(reasons, "primary_nonconvergence")
  if (max(abs(comp$score)) / n >= 1e-6) reasons <- c(reasons, "primary_score_tolerance")
  if (any(!is.finite(comp$probabilities)) || any(comp$probabilities < 0 | comp$probabilities > 1) ||
      max(abs(rowSums(comp$probabilities) - 1)) > 1e-8) reasons <- c(reasons, "invalid_probabilities")
  info <- .matrix_check(comp$H)
  if (!info$valid) reasons <- c(reasons, "information_not_positive_definite")
  if (is.null(alternate$value)) reasons <- c(reasons, "alternate_fit_failed") else {
    alt <- alternate$value; alt_theta <- setNames(.to_natural(alt$par, p), names(theta))
    ac <- ordinal_components(alt_theta, z$X, z$y, K)
    attempts$alternate$convergence <- alt$convergence; attempts$alternate$theta <- alt_theta
    attempts$alternate$logLik <- ac$logLik
    attempts$alternate$max_score_per_n <- if (is.null(ac$score)) Inf else max(abs(ac$score)) / n
    if (alt$convergence != 0 || !is.finite(ac$logLik) || is.null(ac$score) ||
        attempts$alternate$max_score_per_n >= 1e-6) reasons <- c(reasons, "alternate_not_accepted")
    staff <- grep("^(B_c|A_c)", names(theta), value = TRUE)
    if (length(staff) && max(abs(theta[staff] - alt_theta[staff])) > 1e-4)
      reasons <- c(reasons, "staffing_alternate_disagreement")
    if (!is.finite(ac$logLik) || abs(comp$logLik - ac$logLik) > 1e-6 * (1 + abs(comp$logLik)))
      reasons <- c(reasons, "likelihood_alternate_disagreement")
    if (any(!is.finite(alt_theta)) || max(abs(theta - alt_theta)) > 0.01)
      investigation <- c(investigation, "coefficient_growth_or_alternate_instability")
    if (!is.null(ac$H) && .matrix_check(ac$H)$valid && info$valid) {
      se_ratio <- sqrt(diag(chol2inv(chol(ac$H))) / diag(chol2inv(chol(comp$H))))
      if (any(!is.finite(se_ratio)) || max(se_ratio) > 1.1)
        investigation <- c(investigation, "standard_error_growth_under_tighter_tolerance")
      attempts$alternate$information_se_ratio <- se_ratio
    }
  }
  if (min(diff(theta[p + seq_len(K - 1L)])) < 1e-6)
    investigation <- c(investigation, "threshold_gap_below_1e-6")
  if (max(abs(theta)) > 10) investigation <- c(investigation, "absolute_parameter_above_10")
  if (info$condition > 1e10) investigation <- c(investigation, "information_condition_above_1e10")
  colnames(comp$probabilities) <- levels(z$data[[outcome]])
  .finish_fit(z, theta, comp, model_id, attempts, reasons, investigation, "ordinal", alt_theta)
}

binary_components <- function(theta, X, y) {
  eta <- drop(X %*% theta); mu <- plogis(eta)
  scores <- X * as.vector(y - mu)
  H <- crossprod(X, X * as.vector(mu * (1 - mu)))
  ll_i <- ifelse(y == 1, .log_expit(eta), .log_expit(-eta))
  list(logLik = sum(ll_i), logLik_i = ll_i, probabilities = mu,
       scores = scores, score = colSums(scores), H = H)
}

fit_binary <- function(data, outcome, predictors, site = "Site", model_id = "binary") {
  prepared <- .capture(.prepare_model(data, outcome, predictors, site, binary = TRUE))
  if (is.null(prepared$value)) return(.failed_fit(model_id, prepared$error, nrow(data)))
  z <- prepared$value; n <- nrow(z$X)
  primary <- .capture(glm.fit(z$X, z$y, family = binomial(),
    control = glm.control(epsilon = 1e-10, maxit = 500L)))
  attempts <- list(primary = list(solver = "glm.fit IRLS", warnings = primary$warnings, error = primary$error))
  if (is.null(primary$value)) return(.failed_fit(model_id, "primary_fit_failed", n, list(attempts = attempts)))
  m <- primary$value; theta <- setNames(m$coefficients, colnames(z$X))
  if (any(!is.finite(theta))) return(.failed_fit(model_id, "nonfinite_parameters", n, list(attempts = attempts)))
  comp <- binary_components(theta, z$X, z$y)
  alternate <- .capture(optim(theta + 0.03 * sin(seq_along(theta)),
    function(b) -binary_components(b, z$X, z$y)$logLik,
    function(b) -binary_components(b, z$X, z$y)$score,
    method = "BFGS", control = list(maxit = 10000L, reltol = 1e-12)))
  attempts$primary$theta <- theta; attempts$primary$convergence <- m$converged
  attempts$primary$logLik <- comp$logLik
  attempts$alternate <- list(solver = "BFGS tighter tolerance", warnings = alternate$warnings, error = alternate$error)
  reasons <- character(); investigation <- character(); alt_theta <- NULL
  if (!m$converged) reasons <- c(reasons, "primary_nonconvergence")
  if (!is.finite(comp$logLik) || max(abs(comp$score)) / n >= 1e-6)
    reasons <- c(reasons, "nonfinite_likelihood_or_score_tolerance")
  info <- .matrix_check(comp$H)
  if (!info$valid) reasons <- c(reasons, "information_not_positive_definite")
  if (is.null(alternate$value)) reasons <- c(reasons, "alternate_fit_failed") else {
    alt <- alternate$value; alt_theta <- setNames(alt$par, names(theta))
    ac <- binary_components(alt_theta, z$X, z$y)
    attempts$alternate$theta <- alt_theta; attempts$alternate$convergence <- alt$convergence
    attempts$alternate$logLik <- ac$logLik
    attempts$alternate$max_score_per_n <- max(abs(ac$score)) / n
    if (alt$convergence != 0 || !is.finite(ac$logLik) || max(abs(ac$score)) / n >= 1e-6)
      reasons <- c(reasons, "alternate_not_accepted")
    staff <- grep("^(B_c|A_c)", names(theta), value = TRUE)
    if (length(staff) && max(abs(theta[staff] - alt_theta[staff])) > 1e-4)
      reasons <- c(reasons, "staffing_alternate_disagreement")
    if (!is.finite(ac$logLik) || abs(comp$logLik - ac$logLik) > 1e-6 * (1 + abs(comp$logLik)))
      reasons <- c(reasons, "likelihood_alternate_disagreement")
    if (any(!is.finite(alt_theta)) || max(abs(theta - alt_theta)) > 0.01)
      investigation <- c(investigation, "coefficient_growth_or_alternate_instability")
    if (.matrix_check(ac$H)$valid && info$valid) {
      se_ratio <- sqrt(diag(chol2inv(chol(ac$H))) / diag(chol2inv(chol(comp$H))))
      attempts$alternate$information_se_ratio <- se_ratio
      if (any(!is.finite(se_ratio)) || max(se_ratio) > 1.1)
        investigation <- c(investigation, "standard_error_growth_under_tighter_tolerance")
    }
  }
  if (max(abs(theta)) > 10) investigation <- c(investigation, "absolute_parameter_above_10")
  if (info$condition > 1e10) investigation <- c(investigation, "information_condition_above_1e10")
  if (isTRUE(m$boundary) || any(comp$probabilities <= .Machine$double.eps | comp$probabilities >= 1 - .Machine$double.eps))
    investigation <- c(investigation, "boundary_or_separation_warning")
  .finish_fit(z, theta, comp, model_id, attempts, reasons, investigation, "binary", alt_theta)
}

.summary_row <- function(fit, label) {
  data.frame(model_id = fit$model_id, term = label, estimate = NA_real_, se = NA_real_,
    lower = NA_real_, upper = NA_real_, exp_estimate = NA_real_, exp_lower = NA_real_,
    exp_upper = NA_real_, statistic = NA_real_, df1 = NA_real_, df2 = fit$G - 1,
    p_reference = NA_real_, valid = FALSE, status = fit$status, stringsAsFactors = FALSE)
}
.contrast_matrix <- function(fit, weights) {
  if (is.null(names(weights)) || anyDuplicated(names(weights)) ||
      any(!names(weights) %in% names(fit$coef)) || any(!is.finite(weights)))
    stop("Contrast weights must have unique, existing parameter names and finite values.")
  C <- matrix(0, 1L, length(fit$coef), dimnames = list(NULL, names(fit$coef)))
  C[1L, names(weights)] <- weights; C
}

contrast_summary <- function(fit, weights, label = NULL) {
  if (is.null(label)) label <- paste(paste(names(weights), weights, sep = "*"), collapse = "+")
  out <- .summary_row(fit, label); out$df1 <- 1
  if (!length(fit$coef)) return(out)
  C <- .contrast_matrix(fit, weights); out$estimate <- drop(C %*% fit$coef)
  is_threshold <- any(grepl("^threshold:", names(weights)[weights != 0]))
  if (!is_threshold) out$exp_estimate <- exp(out$estimate)
  if (!isTRUE(fit$inference_valid)) {
    out$status <- fit$diagnostics$inference_status %||% fit$status; return(out)
  }
  CV <- C %*% fit$V %*% t(C)
  if (!.matrix_check(CV, 1e10)$valid) {
    out$status <- "contrast_covariance_not_positive_definite_or_ill_conditioned"; return(out)
  }
  out$se <- sqrt(CV[1L, 1L]); out$statistic <- out$estimate / out$se
  margin <- qt(0.975, fit$G - 1L) * out$se
  out$lower <- out$estimate - margin; out$upper <- out$estimate + margin
  if (!is_threshold) { out$exp_lower <- exp(out$lower); out$exp_upper <- exp(out$upper) }
  out$p_reference <- 2 * pt(-abs(out$statistic), fit$G - 1L)
  out$valid <- TRUE; out$status <- "approximate_site_t_reference"; out
}

joint_summary <- function(fit, terms, label = NULL) {
  if (is.null(label)) label <- if (is.character(terms)) paste(terms, collapse = "+") else "joint_contrast"
  out <- .summary_row(fit, label)
  if (!length(fit$coef)) return(out)
  if (is.character(terms)) {
    if (anyDuplicated(terms) || !all(terms %in% names(fit$coef))) stop("Invalid joint contrast terms.")
    C <- matrix(0, length(terms), length(fit$coef), dimnames = list(terms, names(fit$coef)))
    C[cbind(seq_along(terms), match(terms, names(fit$coef)))] <- 1
  } else {
    C <- as.matrix(terms)
    if (is.null(colnames(C)) || anyDuplicated(colnames(C)) || !setequal(colnames(C), names(fit$coef)))
      stop("Joint contrast matrix columns must name every parameter exactly once.")
    C <- C[, names(fit$coef), drop = FALSE]
  }
  h <- nrow(C); out$df1 <- h
  if (h < 1L || h > 2L || any(!is.finite(C)) || qr(t(C))$rank != h) {
    out$status <- "contrast_rows_not_independent_or_dimension_not_prespecified"; return(out)
  }
  if (!isTRUE(fit$inference_valid)) {
    out$status <- fit$diagnostics$inference_status %||% fit$status; return(out)
  }
  CV <- C %*% fit$V %*% t(C)
  if (!.matrix_check(CV, 1e10)$valid) {
    out$status <- "contrast_covariance_not_positive_definite_or_ill_conditioned"; return(out)
  }
  b <- drop(C %*% fit$coef)
  W <- drop(crossprod(b, solve(CV, b)))
  out$statistic <- W / h; out$p_reference <- pf(W / h, h, fit$G - 1L, lower.tail = FALSE)
  out$valid <- TRUE; out$status <- "approximate_site_F_reference"; out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

secondary_holm <- function(p, valid = is.finite(p)) {
  if (length(p) != 2L || length(valid) != 2L) stop("Secondary family contains exactly RQ2 and RQ3.")
  ans <- rep(NA_real_, 2L); names(ans) <- names(p)
  if (isTRUE(all(valid)) && all(is.finite(p)) && all(p >= 0 & p <= 1)) ans <- p.adjust(p, method = "holm")
  ans
}

same_row_comparison <- function(full, reduced) {
  same <- identical(full$rows, reduced$rows) && identical(full$outcome, reduced$outcome) &&
    identical(full$type, reduced$type)
  nested <- same && is.matrix(full$design) && is.matrix(reduced$design) &&
    all(colnames(reduced$design) %in% colnames(full$design)) &&
    isTRUE(all.equal(reduced$design, full$design[, colnames(reduced$design), drop = FALSE], check.attributes = TRUE))
  valid <- isTRUE(full$valid) && isTRUE(reduced$valid) && nested
  delta <- if (valid) full$logLik - reduced$logLik else NA_real_
  if (valid && delta < -1e-6 * (1 + abs(full$logLik))) valid <- FALSE
  data.frame(full_model = full$model_id, reduced_model = reduced$model_id,
    n = full$n, same_rows_outcome = same, nested = nested, valid = valid,
    delta_logLik = delta, apparent_mean_log_score_gain = delta / full$n,
    status = if (valid) "apparent_fit_only" else "invalid_or_unmatched_model_pair")
}

cluster_jackknife_staffing <- function(full, deletedFits, outcome) {
  # Additional approximate uncertainty sensitivity. Retain the full-fit
  # staffing estimates as the center. Nuisance parameters change with the
  # deleted site and are deliberately excluded from this covariance.
  if (!outcome %in% c("Q6", "Q4")) stop("outcome must be Q6 or Q4.")
  staff <- c("B_c", "A_c")
  C <- if (outcome == "Q6") rbind(B_c = c(1, 0), A_c = c(0, 1), A_minus_B = c(-1, 1)) else diag(2)
  if (outcome == "Q4") rownames(C) <- staff
  colnames(C) <- staff
  n_valid <- sum(vapply(deletedFits, function(f) isTRUE(f$valid), logical(1)))
  out <- data.frame(model_id = full$model_id, outcome = outcome, term = rownames(C),
    estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_,
    exp_estimate = NA_real_, exp_lower = NA_real_, exp_upper = NA_real_,
    df = 7L, valid = FALSE, status = "jackknife_unavailable",
    n_deletions = length(deletedFits), n_valid_deletions = n_valid,
    stringsAsFactors = FALSE)
  Vjk <- matrix(NA_real_, 2L, 2L, dimnames = list(staff, staff))
  reasons <- character()
  if (!isTRUE(full$valid)) reasons <- c(reasons, "invalid_full_fit")
  if (!all(staff %in% names(full$coef)) || any(!is.finite(full$coef[staff])))
    reasons <- c(reasons, "full_staffing_coefficients_unavailable")
  if (isTRUE(full$valid) && !"full_staffing_coefficients_unavailable" %in% reasons) {
    out$estimate <- drop(C %*% full$coef[staff]); out$exp_estimate <- exp(out$estimate)
  }
  full_sites <- sort(unique(as.character(full$site)))
  if (!identical(full_sites, as.character(1:8))) reasons <- c(reasons, "full_fit_requires_sites_1_to_8")
  if (length(full$rows) != length(full$site) || anyDuplicated(full$rows))
    reasons <- c(reasons, "full_fit_row_identifiers_invalid")
  if (length(deletedFits) != 8L) reasons <- c(reasons, "exactly_eight_deletions_required")
  if (n_valid != 8L) reasons <- c(reasons, "one_or_more_deletion_fits_invalid")
  omitted <- rep(NA_character_, length(deletedFits))
  deletion_status <- vapply(deletedFits, function(f) f$status %||% "status_missing", character(1))
  for (j in seq_along(deletedFits)) {
    f <- deletedFits[[j]]
    kept <- sort(unique(as.character(f$site)))
    gone <- setdiff(full_sites, kept)
    if (length(gone) != 1L || length(kept) != 7L || !all(kept %in% full_sites)) {
      reasons <- c(reasons, "deletion_site_membership_invalid"); next
    }
    omitted[j] <- gone
    wanted <- full$rows[as.character(full$site) != gone]
    if (anyDuplicated(f$rows) || !setequal(wanted, f$rows) || length(wanted) != length(f$rows))
      reasons <- c(reasons, "deletion_rows_do_not_match_full_fit")
    else if (!identical(f$outcome, full$outcome[match(f$rows, full$rows)]))
      reasons <- c(reasons, "deletion_outcome_does_not_match_full_fit")
    if (!identical(f$type, full$type)) reasons <- c(reasons, "deletion_model_type_differs")
    if (!all(staff %in% names(f$coef)) || any(!is.finite(f$coef[staff])))
      reasons <- c(reasons, "deletion_staffing_coefficients_unavailable")
  }
  if (!identical(sort(omitted), as.character(1:8))) reasons <- c(reasons, "each_site_must_be_deleted_once")
  if (!length(reasons)) {
    estimates <- do.call(rbind, lapply(deletedFits, function(f) unname(f$coef[staff])))
    delta <- sweep(estimates, 2L, unname(full$coef[staff]), "-")
    Vjk <- (7 / 8) * crossprod(delta); dimnames(Vjk) <- list(staff, staff)
    for (j in seq_len(nrow(C))) {
      variance <- C[j, , drop = FALSE] %*% Vjk %*% t(C[j, , drop = FALSE])
      if (!.matrix_check(variance, 1e10)$valid) {
        out$status[j] <- "jackknife_contrast_variance_not_positive_finite"; next
      }
      out$se[j] <- sqrt(variance[1L, 1L]); margin <- qt(0.975, 7) * out$se[j]
      out$lower[j] <- out$estimate[j] - margin; out$upper[j] <- out$estimate[j] + margin
      out$exp_lower[j] <- exp(out$lower[j]); out$exp_upper[j] <- exp(out$upper[j])
      out$valid[j] <- TRUE
      out$status[j] <- "approximate_full_centered_cluster_jackknife_t7_sensitivity"
    }
  } else out$status <- paste(unique(reasons), collapse = ";")
  attr(out, "Vjk") <- Vjk
  attr(out, "deletion_status") <- deletion_status
  attr(out, "omitted_site") <- omitted
  out
}
