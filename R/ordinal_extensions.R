# Specification 1.4, priority 2. Sourcing defines functions and runs nothing.
# Entry: run_ordinal_extensions(data, source_run, outdir).
# data is the FULL decoded frame returned by read_study_data(), not a CC frame.
# Source data_prep.R and model_engine.R first. source_run contains the accepted
# original fits. outdir must not already exist. All new results are point-only.

.oe_capture <- function(expr) .capture(expr)
.oe_common <- c("L_c", "M_c", "B_c", "A_c")
.oe_numeric <- c("L_c", "E_c", "R_c", "M_c", "B_c", "A_c")
.oe_factor_predictors <- c("L_c", "E_c", "R_c", "M_c", "B_f", "A_f")
.oe_eps <- 1e-8

.oe_data <- function(data) {
  d <- data
  d$Site <- factor(as.character(d$Site), levels = as.character(1:8))
  if (anyNA(d$Site)) stop("Unsupported site in ordinal extension.")
  contrasts(d$Site) <- contr.treatment(levels(d$Site), base = 1L)
  d
}

.oe_factor_data <- function(data) {
  d <- .oe_data(data)
  for (v in c("B", "A")) {
    if (anyNA(d[[v]]) || any(!d[[v]] %in% 1:5)) stop("Invalid staffing category.")
    d[[paste0(v, "_f")]] <- factor(d[[v]], levels = c(3, 1, 2, 4, 5))
    contrasts(d[[paste0(v, "_f")]]) <- contr.treatment(c("3", "1", "2", "4", "5"), base = 1L)
  }
  d
}

# Natural order is tau[1:3], gamma_E[1:3], gamma_R[1:3], common slopes.
oe_ppo_design <- function(data, check_rank = TRUE) {
  d <- .oe_data(data)
  required <- c("E_c", "R_c", .oe_common, "Site", "Y")
  if (!all(required %in% names(d)) || anyNA(d[, required])) stop("Incomplete PPO design.")
  y <- as.integer(ordered(d$Y, levels = 1:4))
  if (anyNA(y) || (check_rank && any(tabulate(y, nbins = 4L) == 0))) stop("Invalid four-category outcome.")
  X <- model.matrix(reformulate(c("Site", .oe_common)), d)[, -1L, drop = FALSE]
  q <- 9L + ncol(X)
  nm <- c(paste0("tau", 1:3), paste0("E", 1:3), paste0("R", 1:3), colnames(X))
  W <- lapply(1:3, function(k) {
    w <- matrix(0, nrow(d), q, dimnames = list(NULL, nm))
    w[, k] <- 1; w[, 3L + k] <- -d$E_c; w[, 6L + k] <- -d$R_c
    w[, 9L + seq_len(ncol(X))] <- -X
    w
  })
  pairs <- unique(rbind(data.frame(E_c = d$E_c, R_c = d$R_c), data.frame(E_c = 0, R_c = 0)))
  C <- do.call(rbind, lapply(1:2, function(k) {
    v <- matrix(0, nrow(pairs), q, dimnames = list(NULL, nm))
    v[, k + 1L] <- 1; v[, k] <- -1
    v[, 3L + k + 1L] <- -pairs$E_c; v[, 3L + k] <- pairs$E_c
    v[, 6L + k + 1L] <- -pairs$R_c; v[, 6L + k] <- pairs$R_c
    v
  }))
  if (check_rank && (q != 20L || qr(do.call(rbind, W))$rank != q)) stop("PPO stacked design is not rank 20.")
  list(data = d, X = X, W = W, C = C, pairs = pairs, y = y, q = q,
       names = nm, n = nrow(d), site = d$Site, rows = as.character(d$Observation))
}

oe_ppo_components <- function(theta, z, derivatives = TRUE) {
  bad <- function() list(logLik = -Inf, probabilities = matrix(NA_real_, z$n, 4L))
  if (length(theta) != z$q || any(!is.finite(theta))) return(bad())
  a <- vapply(z$W, function(w) drop(w %*% theta), numeric(z$n))
  gaps <- a[, 2:3, drop = FALSE] - a[, 1:2, drop = FALSE]
  # Every category at every observed covariate row must be valid, including
  # categories that were not observed for that individual.
  if (any(!is.finite(a)) || any(gaps <= 0) || any(drop(z$C %*% theta) <= 0)) return(bad())
  lp <- cbind(.log_expit(a[, 1]),
    .log_expit(a[, 2]) + .log_expit(-a[, 1]) + .log_one_minus_exp(-gaps[, 1]),
    .log_expit(a[, 3]) + .log_expit(-a[, 2]) + .log_one_minus_exp(-gaps[, 2]),
    .log_expit(-a[, 3]))
  ll <- lp[cbind(seq_len(z$n), z$y)]
  ans <- list(logLik = sum(ll), logLik_i = ll, probabilities = exp(lp),
              cumulative_predictors = a, gaps = gaps)
  colnames(ans$probabilities) <- as.character(1:4)
  if (!derivatives) return(ans)
  wa <- wb <- matrix(0, z$n, z$q, dimnames = list(NULL, z$names))
  Fa <- numeric(z$n); Fb <- rep(1, z$n); r <- numeric(z$n)
  for (k in 1:3) {
    lo <- which(z$y == k + 1L); hi <- which(z$y == k)
    if (length(lo)) { wa[lo, ] <- z$W[[k]][lo, , drop = FALSE]; Fa[lo] <- plogis(a[lo, k]) }
    if (length(hi)) { wb[hi, ] <- z$W[[k]][hi, , drop = FALSE]; Fb[hi] <- plogis(a[hi, k]) }
  }
  mid <- which(z$y %in% 2:3)
  if (length(mid)) r[mid] <- 1 / expm1(a[cbind(mid, z$y[mid])] - a[cbind(mid, z$y[mid] - 1L)])
  scores <- wa * (-Fa - r) + wb * (1 - Fb + r)
  wd <- wb - wa
  H <- crossprod(wa, wa * (Fa * (1 - Fa))) + crossprod(wb, wb * (Fb * (1 - Fb))) +
    crossprod(wd, wd * (r * (1 + r)))
  dimnames(H) <- list(z$names, z$names)
  ans$scores <- scores; ans$score <- colSums(scores); ans$H <- (H + t(H)) / 2
  ans
}

oe_po_to_ppo <- function(fit, z) {
  if (!isTRUE(fit$valid) || length(fit$coef) != 16L) stop("Accepted original PO initializer required.")
  J <- matrix(0, z$q, length(fit$coef), dimnames = list(z$names, names(fit$coef)))
  th <- grep("^threshold:", names(fit$coef))
  if (length(th) != 3L || !all(c("E_c", "R_c", colnames(z$X)) %in% names(fit$coef))) stop("Unexpected original PO parameter names.")
  J[cbind(1:3, th)] <- 1
  J[4:6, "E_c"] <- 1; J[7:9, "R_c"] <- 1
  J[cbind(9L + seq_len(ncol(z$X)), match(colnames(z$X), names(fit$coef)))] <- 1
  list(theta = setNames(drop(J %*% fit$coef), z$names), J = J)
}

.oe_hull_reference <- function(e, r, tolerance = 1e-10) {
  u <- unique(cbind(e, r)); target <- c(0, 0)
  if (any(rowSums(abs(u)) <= tolerance)) return(TRUE)
  if (nrow(u) == 1L) return(FALSE)
  if (qr(sweep(u, 2, u[1, ]))$rank < 2L) {
    j <- which.max(rowSums((u - matrix(u[1, ], nrow(u), 2L, byrow = TRUE))^2))
    v <- u[j, ] - u[1, ]; den <- sum(v^2)
    if (den == 0) return(FALSE)
    t <- sum((target - u[1, ]) * v) / den
    coordinate <- drop(sweep(u, 2, u[1, ]) %*% v) / den
    return(t >= min(coordinate) - tolerance && t <= max(coordinate) + tolerance &&
             max(abs(u[1, ] + t * v)) <= tolerance)
  }
  h <- u[chull(u), , drop = FALSE]; nxt <- h[c(2:nrow(h), 1L), , drop = FALSE]
  cross <- (nxt[, 1] - h[, 1]) * (-h[, 2]) - (nxt[, 2] - h[, 2]) * (-h[, 1])
  all(cross >= -tolerance) || all(cross <= tolerance)
}

# Mapping is checked against actual package X/NOM design matrices. The public
# clm control explicitly uses positive nominal and negative location signs.
.oe_clm_map <- function(model, z) {
  dm <- model.matrix(model, type = "design")
  if (!identical(colnames(dm$X), c("(Intercept)", colnames(z$X))) ||
      !identical(colnames(dm$NOM), c("(Intercept)", "E_c", "R_c"))) stop("Unexpected clm design ordering.")
  if (max(abs(dm$X[, -1L, drop = FALSE] - z$X)) > 1e-12 ||
      max(abs(dm$NOM[, -1L, drop = FALSE] - cbind(z$data$E_c, z$data$R_c))) > 1e-12)
    stop("clm design does not match natural design.")
  phi <- coef(model)
  if (length(phi) != z$q || any(!is.finite(phi)) || length(model$alpha) != 9L || length(model$beta) != 11L)
    stop("Dropped or unexpected clm coefficient.")
  signs <- c(rep(1, 3L), rep(-1, 6L), rep(1, 11L))
  P <- diag(signs); dimnames(P) <- list(z$names, names(phi))
  theta <- setNames(drop(P %*% phi), z$names)
  for (k in 1:3) {
    v <- matrix(0, z$n, z$q)
    for (j in 1:3) v[, (j - 1L) * 3L + k] <- dm$NOM[, j]
    v[, 10:20] <- -dm$X[, -1L, drop = FALSE]
    if (max(abs(v %*% phi - z$W[[k]] %*% theta)) > 1e-10) stop("Nominal sign mapping failed.")
  }
  list(theta = theta, P = P, package_names = names(phi), design = dm)
}

.oe_assess <- function(theta, z) {
  c <- oe_ppo_components(theta, z)
  if (!is.finite(c$logLik) || is.null(c$H)) return(list(ok = FALSE, reasons = "nonfinite_likelihood_or_derivatives", comp = c))
  info <- .matrix_check(c$H)
  slack <- drop(z$C %*% theta) - .oe_eps
  reasons <- character()
  if (max(abs(c$score)) / z$n >= 1e-6) reasons <- c(reasons, "unrestricted_score_tolerance")
  if (!info$valid || info$condition >= 1e10) reasons <- c(reasons, "information_not_accepted")
  if (any(!is.finite(c$probabilities)) || any(c$probabilities <= 0 | c$probabilities > 1) ||
      max(abs(rowSums(c$probabilities) - 1)) > 1e-8) reasons <- c(reasons, "invalid_probabilities")
  if (any(!is.finite(c$scores)) || any(!is.finite(c$H))) reasons <- c(reasons, "nonfinite_derivatives")
  if (min(slack) <= 1e-6) reasons <- c(reasons, "boundary_or_near_boundary_unavailable")
  if (max(abs(theta)) > 10) reasons <- c(reasons, "absolute_parameter_above_10_investigation")
  list(ok = !length(reasons), reasons = reasons, comp = c, information = info,
       min_slack = min(slack), active_constraints = which(slack <= 1e-6),
       max_score_per_n = max(abs(c$score)) / z$n)
}

.oe_compare <- function(theta, candidate, z, label) {
  if (is.null(candidate) || any(!is.finite(candidate))) return(paste0(label, "_unavailable"))
  a <- oe_ppo_components(theta, z); b <- oe_ppo_components(candidate, z)
  if (!is.finite(a$logLik) || !is.finite(b$logLik)) return(paste0(label, "_nonfinite"))
  reasons <- character()
  if (max(abs(theta[c("B_c", "A_c")] - candidate[c("B_c", "A_c")])) > 1e-4) reasons <- c(reasons, "staffing")
  if (max(abs(theta - candidate)) > 0.01) reasons <- c(reasons, "all_parameters")
  if (abs(a$logLik - b$logLik) > 1e-6 * (1 + abs(a$logLik))) reasons <- c(reasons, "likelihood")
  if (max(abs(a$probabilities - b$probabilities)) > 1e-5) reasons <- c(reasons, "probabilities")
  # Checking all original covariate rows under each fixed staffing change
  # covers every supported display subset without selecting a fitted result.
  for (item in c("B", "A")) for (value in c(3, 4)) {
    d <- z$data; d[[item]] <- value; d[[paste0(item, "_c")]] <- value - 3
    zz <- oe_ppo_design(d, FALSE)
    ap <- oe_ppo_components(theta, zz, FALSE)$probabilities
    bp <- oe_ppo_components(candidate, zz, FALSE)$probabilities
    if (any(!is.finite(ap)) || any(!is.finite(bp)) || max(abs(ap - bp)) > 1e-5)
      reasons <- c(reasons, "scenario_probabilities")
  }
  reasons <- unique(reasons)
  paste0(label, "_", reasons, "_disagreement")[seq_along(reasons)]
}

.oe_constrained_attempt <- function(start, z, mu, reltol) {
  if (any(drop(z$C %*% start) <= .oe_eps)) return(list(theta = NULL, error = "infeasible_start"))
  f <- function(p) {
    if (any(!is.finite(p)) || any(drop(z$C %*% p) <= .oe_eps)) return(Inf)
    -oe_ppo_components(p, z, FALSE)$logLik
  }
  gr <- function(p) -oe_ppo_components(p, z)$score
  res <- .oe_capture(constrOptim(start, f, gr, ui = z$C, ci = rep(.oe_eps, nrow(z$C)),
    mu = mu, method = "BFGS", control = list(maxit = 10000L, reltol = reltol),
    outer.iterations = 100L, outer.eps = 1e-10))
  record <- list(start = start, mu = mu, reltol = reltol, result = res,
                 theta = if (is.null(res$value)) NULL else setNames(res$value$par, z$names))
  if (!is.null(record$theta)) {
    record$assessment <- .oe_assess(record$theta, z)
    # Removing the numerical barrier is allowed only inside the identical
    # feasible set. The original barrier solve is retained in this record.
    pol <- .oe_capture(optim(record$theta, f, gr, method = "BFGS",
                            control = list(maxit = 10000L, reltol = 1e-13)))
    record$polish <- pol
    if (!is.null(pol$value) && isTRUE(pol$value$convergence == 0L)) {
      tp <- setNames(pol$value$par, z$names); cp <- oe_ppo_components(tp, z)
      old <- record$assessment$comp
      if (is.finite(cp$logLik) && cp$logLik >= old$logLik - 100 * .Machine$double.eps * (1 + abs(old$logLik)) &&
          !is.null(cp$score) && max(abs(cp$score)) <= max(abs(old$score))) {
        record$theta <- tp; record$assessment <- .oe_assess(tp, z); record$selected <- "interior_polish"
      }
    }
    if (is.null(record$selected)) record$selected <- "barrier"
    record$convergence <- if (record$selected == "interior_polish") pol$value$convergence else res$value$convergence
  }
  record
}

fit_ppo_extension <- function(data, original_fit, model_id = "ppo_ER") {
  z <- oe_ppo_design(data)
  base <- oe_po_to_ppo(original_fit, z)$theta
  if (any(drop(z$C %*% base) <= .oe_eps)) stop("Original fit is not a strictly feasible PPO initializer.")
  ctl <- ordinal::clm.control(sign.location = "negative", sign.nominal = "positive",
    maxIter = 1000L, gradTol = 1e-7, relTol = 1e-10, maxLineIter = 50L)
  # package alpha adds E/R, whereas the natural coefficients subtract them.
  package_start <- unname(base * c(rep(1, 3), rep(-1, 6), rep(1, 11)))
  primary <- .oe_capture(ordinal::clm(Y ~ Site + L_c + M_c + B_c + A_c,
    nominal = ~ E_c + R_c, data = z$data, link = "logit", threshold = "flexible",
    start = package_start, model = TRUE, control = ctl))
  mapping <- if (is.null(primary$value)) list(value = NULL, error = primary$error) else
    .oe_capture(.oe_clm_map(primary$value, z))
  attempts <- list(package = primary, mapping = mapping)
  reasons <- character(); theta <- NULL; assess <- NULL
  if (is.null(mapping$value)) reasons <- c(reasons, "package_fit_or_mapping_failed") else {
    theta <- mapping$value$theta; assess <- .oe_assess(theta, z)
    reasons <- c(reasons, assess$reasons)
    code <- primary$value$convergence
    if (is.list(code)) code <- code$code
    if (!isTRUE(code == 0L)) reasons <- c(reasons, "package_nonconvergence")
    if (is.finite(assess$comp$logLik)) {
      if (abs(as.numeric(logLik(primary$value)) - assess$comp$logLik) > 1e-7 * (1 + abs(assess$comp$logLik)))
        reasons <- c(reasons, "package_likelihood_mismatch")
      P <- mapping$value$P
      Hp <- P %*% primary$value$Hessian %*% t(P)
      gp <- -drop(P %*% primary$value$gradient)
      attempts$package_derivative_error <- c(score = max(abs(gp - assess$comp$score)) / (1 + max(abs(gp))),
        information = max(abs(Hp - assess$comp$H)) / (1 + max(abs(Hp))))
      if (any(!is.finite(attempts$package_derivative_error)) ||
          attempts$package_derivative_error["score"] > 1e-5 || attempts$package_derivative_error["information"] > 1e-4)
        reasons <- c(reasons, "package_derivative_mismatch")
      pd <- z$data; pd$Y <- NULL
      pp <- .oe_capture(predict(primary$value, newdata = pd, type = "prob")$fit)
      attempts$package_all_category_predictions <- pp
      if (is.null(pp$value) || !identical(dim(pp$value), c(z$n, 4L)) ||
          any(!is.finite(pp$value)) || max(abs(pp$value - assess$comp$probabilities)) > 1e-5)
        reasons <- c(reasons, "package_all_category_prediction_mismatch")
    }
  }
  delta <- 0.03 * sin(seq_along(base)); halvings <- 0L
  while (any(drop(z$C %*% (base + delta)) <= .oe_eps) && halvings < 60L) {
    delta <- delta / 2; halvings <- halvings + 1L
  }
  attempts$perturbation_halvings <- halvings
  starts <- list(original = base, perturbed = base + delta)
  for (nm in names(starts)) {
    ca <- .oe_constrained_attempt(starts[[nm]], z, 1e-4, 1e-12)
    attempts[[paste0("constrained_", nm)]] <- ca
    selected <- ca
    # The trigger is evaluated before any permitted score polish as well as
    # after it, so barrier dependence near a boundary is never hidden.
    rawtheta <- if (!is.null(ca$result$value)) setNames(ca$result$value$par, z$names) else NULL
    rawcheck <- if (is.null(rawtheta)) NULL else .oe_assess(rawtheta, z)
    retry <- !is.null(rawcheck) && is.finite(rawcheck$comp$logLik) &&
      (rawcheck$max_score_per_n >= 1e-6 || rawcheck$min_slack < 1e-4)
    if (retry) {
      tight <- .oe_constrained_attempt(if (!is.null(ca$theta)) ca$theta else starts[[nm]], z, 1e-6, 1e-13)
      attempts[[paste0("constrained_", nm, "_tight")]] <- tight
      if (!is.null(ca$theta) && !is.null(tight$theta)) {
        reasons <- c(reasons, .oe_compare(ca$theta, tight$theta, z, paste0(nm, "_barrier")))
        h1 <- ca$assessment$comp$H; h2 <- tight$assessment$comp$H
        if (!is.null(h1) && !is.null(h2) && .matrix_check(h1)$valid && .matrix_check(h2)$valid) {
          ratio <- sqrt(diag(chol2inv(chol(h2))) / diag(chol2inv(chol(h1))))
          attempts[[paste0(nm, "_information_se_ratio")]] <- ratio
          if (any(!is.finite(ratio)) || max(ratio) > 1.1) reasons <- c(reasons, "information_se_growth")
        }
      }
      selected <- tight
    }
    attempts[[paste0("selected_", nm)]] <- selected$selected
    if (is.null(selected$theta) || !isTRUE(selected$convergence == 0L) || !isTRUE(selected$assessment$ok))
      reasons <- c(reasons, paste0(nm, "_crosscheck_not_accepted"))
    if (!is.null(theta)) reasons <- c(reasons, .oe_compare(theta, selected$theta, z, nm))
  }
  if (!is.null(assess) && is.finite(assess$comp$logLik) &&
      assess$comp$logLik < original_fit$logLik - 1e-6 * (1 + abs(original_fit$logLik))) reasons <- c(reasons, "nested_likelihood_loss")
  reasons <- unique(reasons)
  valid <- !length(reasons)
  list(model_id = model_id, type = "ppo_ER", valid = valid, inference_valid = FALSE,
    status = if (valid) "accepted; point estimates only" else paste(reasons, collapse = ";"),
    coef = if (valid) theta else numeric(), candidate_theta = theta,
    probabilities = if (valid) assess$comp$probabilities else NULL,
    scores = if (valid) assess$comp$scores else NULL, H = if (valid) assess$comp$H else NULL,
    logLik = if (valid) assess$comp$logLik else NA_real_, n = z$n, G = length(unique(z$site)),
    rows = z$rows, site = z$site, design = z$X, predictors = .oe_common,
    diagnostics = list(reasons = reasons, assessment = assess, attempts = attempts,
      constraints = z$C, support_pairs = z$pairs,
      reference_inside_observed_hull = .oe_hull_reference(z$data$E_c, z$data$R_c),
      reference_gap = if (is.null(theta)) rep(NA_real_, 2L) else diff(theta[1:3]),
      package_version = as.character(utils::packageVersion("ordinal"))))
}

.oe_po_probabilities <- function(fit, data, predictors) {
  d <- if ("B_f" %in% predictors) .oe_factor_data(data) else .oe_data(data)
  X <- model.matrix(reformulate(c("Site", predictors)), d)[, -1L, drop = FALSE]
  if (!identical(colnames(X), colnames(fit$design))) stop("Prediction design differs from fitted PO design.")
  theta <- fit$coef
  ordinal_components(theta, X, as.integer(d$Y), 4L, FALSE)$probabilities
}

oe_extension_predict <- function(fit, data) {
  if (!isTRUE(fit$valid)) stop("Prediction requested from unavailable model.")
  p <- if (identical(fit$type, "ppo_ER")) oe_ppo_components(fit$coef, oe_ppo_design(data, FALSE), FALSE)$probabilities else
    .oe_po_probabilities(fit, data, fit$predictors)
  if (any(!is.finite(p)) || any(p <= 0 | p > 1) || max(abs(rowSums(p) - 1)) > 1e-8)
    stop("Invalid all-category prediction.")
  colnames(p) <- as.character(1:4); p
}

oe_staffing_support <- function(data, outcome) {
  d <- .oe_data(data)
  cells <- as.data.frame(table(Site = factor(d$Site, levels = 1:8), B = factor(d$B, levels = 1:5),
                              A = factor(d$A, levels = 1:5)), stringsAsFactors = FALSE)
  names(cells)[4] <- "n"
  key <- function(g, b, a) paste(g, b, a, sep = ":")
  counts <- setNames(cells$n, key(cells$Site, cells$B, cells$A))
  b3 <- unname(counts[key(d$Site, 3, d$A)]); b4 <- unname(counts[key(d$Site, 4, d$A)])
  a3 <- unname(counts[key(d$Site, d$B, 3)]); a4 <- unname(counts[key(d$Site, d$B, 4)])
  eligible <- list(B = b3 >= 5L & b4 >= 5L, A = a3 >= 5L & a4 >= 5L)
  if (outcome == 6) eligible$AB <- eligible$A & eligible$B
  status <- do.call(rbind, lapply(names(eligible), function(v) {
    x <- eligible[[v]]; n <- sum(x); g <- length(unique(d$Site[x]))
    data.frame(outcome = outcome, population = v, n_original = nrow(d), n = n,
      excluded = nrow(d) - n, sites = g, supported = n >= 30L && g >= 3L,
      status = if (n >= 30L && g >= 3L) "supported by prespecified display rule" else "insufficient prespecified staffing support")
  }))
  rows <- data.frame(outcome = outcome, Observation = as.character(d$Observation), Site = as.character(d$Site),
    B = d$B, A = d$A, n_B3 = b3, n_B4 = b4, n_A3 = a3, n_A4 = a4,
    eligible_B = eligible$B, eligible_A = eligible$A)
  if (outcome == 6) rows$eligible_AB <- eligible$AB
  composition <- do.call(rbind, lapply(names(eligible), function(v) {
    t <- table(factor(d$Site[eligible[[v]]], levels = 1:8))
    data.frame(outcome = outcome, population = v, Site = names(t), n = as.integer(t))
  }))
  cells$outcome <- outcome
  list(eligible = eligible, status = status, rows = rows, composition = composition, cells = cells)
}

oe_probability_summary <- function(fits, data, outcome, support = oe_staffing_support(data, outcome)) {
  probs <- contrasts <- statuses <- list()
  categories <- c("disagreement", "neutral", "agree", "strongly_agree", "agree_or_strongly_agree")
  scenario <- function(fit, d, item, level) {
    d[[item]] <- level; d[[paste0(item, "_c")]] <- level - 3
    p <- oe_extension_predict(fit, d)
    c(colMeans(p), mean(rowSums(p[, 3:4, drop = FALSE])))
  }
  for (id in names(fits)) for (pop in names(support$eligible)) {
    f <- fits[[id]]; sr <- support$status[support$status$population == pop, ]
    status <- if (!isTRUE(f$valid)) paste("model unavailable:", f$status) else if (!sr$supported) sr$status else "accepted; point estimates only"
    result <- NULL
    if (isTRUE(f$valid) && sr$supported) {
      dd <- data[support$eligible[[pop]], , drop = FALSE]
      result <- .oe_capture({
        items <- if (pop == "AB") c("A", "B") else pop
        changes <- list(); distributions <- list()
        for (item in items) {
          p3 <- scenario(f, dd, item, 3); p4 <- scenario(f, dd, item, 4)
          changes[[item]] <- p4 - p3
          distributions[[item]] <- rbind(
            data.frame(item = item, level = 3L, category = categories, probability = p3),
            data.frame(item = item, level = 4L, category = categories, probability = p4))
        }
        cc <- if (pop == "AB") changes$A - changes$B else changes[[pop]]
        if (abs(sum(cc[1:4])) > 1e-8) stop("Category contrasts do not sum to zero.")
        list(distributions = do.call(rbind, distributions), contrast = cc)
      })
      if (is.null(result$value)) status <- paste("prediction unavailable:", result$error)
    }
    statuses[[length(statuses) + 1L]] <- data.frame(model_id = id, outcome = outcome, population = pop,
      n = sr$n, sites = sr$sites, valid = !is.null(result$value), status = status)
    if (!is.null(result$value)) {
      pp <- result$value$distributions; pp$model_id <- id; pp$outcome <- outcome; pp$population <- pop
      pp$n <- sr$n; pp$sites <- sr$sites; probs[[length(probs) + 1L]] <- pp
      contrasts[[length(contrasts) + 1L]] <- data.frame(model_id = id, outcome = outcome, population = pop,
        contrast = if (pop == "AB") "A_minus_B_3_to_4_common_population" else paste0(pop, "_3_to_4"),
        category = categories, difference = result$value$contrast, percentage_points = 100 * result$value$contrast,
        n = sr$n, sites = sr$sites, status = "point estimate only")
    }
  }
  bind <- function(x) if (length(x)) do.call(rbind, x) else data.frame(status = "No supported accepted predictions.")
  list(probabilities = bind(probs), contrasts = bind(contrasts), status = bind(statuses))
}

oe_refuse_joint_factor_wald <- function(...) stop("Eight staffing restrictions exceed the seven-rank site bound; no factor Wald test is specified.")

oe_factor_alternate_reasons <- function(fit) {
  terms <- paste0(rep(c("B_f", "A_f"), each = 4L), rep(c(1, 2, 4, 5), 2L))
  if (!all(terms %in% names(fit$coef)) || !all(terms %in% names(fit$alternate)) ||
      any(!is.finite(fit$coef[terms])) || any(!is.finite(fit$alternate[terms])))
    return("factor_staffing_alternate_unavailable")
  if (max(abs(fit$coef[terms] - fit$alternate[terms])) > 1e-4)
    return("factor_staffing_alternate_disagreement")
  character()
}

.oe_point_results <- function(fit, outcome) {
  terms <- if (identical(fit$type, "five_level_PO")) paste0(rep(c("B_f", "A_f"), each = 4L), rep(c(1, 2, 4, 5), 2L)) else c("B_c", "A_c")
  if (outcome == 6 && !identical(fit$type, "five_level_PO")) terms <- c(terms, "A_minus_B")
  estimates <- vapply(terms, function(t) {
    if (!isTRUE(fit$valid)) return(NA_real_)
    if (t == "A_minus_B") unname(fit$coef["A_c"] - fit$coef["B_c"]) else unname(fit$coef[t])
  }, numeric(1))
  data.frame(model_id = fit$model_id, outcome = outcome, term = terms, estimate = estimates,
    exp_estimate = exp(estimates), n = fit$n, sites = fit$G, valid = fit$valid, status = fit$status,
    stringsAsFactors = FALSE)
}

oe_validate_source_fit <- function(f, z, expected_model_id) {
  lev <- levels(z$data$Y)
  expected_names <- c(colnames(z$X), paste0("threshold:", lev[-length(lev)], "|", lev[-1L]))
  metadata <- f$analysis_metadata
  if (!isTRUE(f$valid) || !identical(f$model_id, expected_model_id) || !identical(f$type, "ordinal") ||
      !isTRUE(f$n == nrow(z$X)) || !isTRUE(f$G == 8L) ||
      !identical(names(f$coef), expected_names) || any(!is.finite(f$coef)) ||
      !is.list(metadata) || !identical(metadata$outcome, "Y") || !identical(metadata$kind, "ordinal") ||
      !identical(metadata$predictors, .oe_numeric) ||
      !identical(as.character(metadata$row_identifiers), z$rows) ||
      !identical(as.character(f$rows), z$rows) ||
      !identical(as.integer(f$outcome), z$y) || !identical(as.character(f$site), as.character(z$site)) ||
      !identical(dim(f$design), dim(z$X)) || !identical(colnames(f$design), colnames(z$X)) ||
      max(abs(f$design - z$X)) > 1e-12) stop("Saved main fit/data provenance mismatch.")
  check <- ordinal_components(f$coef, z$X, z$y, 4L, FALSE)
  if (!is.finite(f$logLik) || !is.finite(check$logLik) ||
      abs(f$logLik - check$logLik) > 1e-7 * (1 + abs(check$logLik)))
    stop("Saved main likelihood does not match its ordered coefficients and design.")
  invisible(TRUE)
}

.oe_source_fit <- function(source_run, outcome, data) {
  expected_id <- paste0("q", outcome, "_main_full")
  path <- file.path(source_run, "fits", paste0(expected_id, ".rds"))
  if (!file.exists(path)) stop("Missing original main fit: ", path)
  f <- readRDS(path)
  z <- .prepare_model(data, "Y", .oe_numeric, "Site")
  oe_validate_source_fit(f, z, expected_id)
  f$predictors <- .oe_numeric
  f$inference_valid <- FALSE # this module never reissues primary inference
  attr(f, "source_path") <- normalizePath(path, winslash = "/")
  f
}

run_ordinal_extensions <- function(data, source_run, outdir) {
  required <- c("read_study_data", "derive_scores", "complete_case_rows", "fit_ordinal", "ordinal_components", ".prepare_model")
  if (!all(vapply(required, exists, logical(1), envir = environment(), mode = "function", inherits = TRUE))) stop("Source data_prep.R and model_engine.R first.")
  if (file.exists(outdir) || dir.exists(outdir)) stop("Ordinal output directory must be new.")
  if (!dir.create(outdir, recursive = TRUE)) stop("Could not create ordinal output directory.")
  dir.create(file.path(outdir, "fits"))
  pkg <- .oe_capture(requireNamespace("ordinal", quietly = TRUE))
  saveRDS(pkg, file.path(outdir, "ordinal_package_load.rds"))
  if (!isTRUE(pkg$value)) stop("ordinal package unavailable; inspect package load record.")
  write.csv(data.frame(package = c("R", "ordinal", "MASS"), version = c(as.character(getRversion()),
    as.character(utils::packageVersion("ordinal")), as.character(utils::packageVersion("MASS")))),
    file.path(outdir, "package_versions.csv"), row.names = FALSE)
  fits_all <- list(); estimates <- statuses <- supports <- predictions <- list(); source_files <- character()
  for (o in c(6, 4)) {
    d <- derive_scores(data, o); d <- d[complete_case_rows(data, o), , drop = FALSE]
    d <- .oe_data(d)
    original <- .oe_source_fit(source_run, o, d)
    source_files <- c(source_files, attr(original, "source_path"))
    ppo_id <- paste0("q", o, "_ppo_ER"); fac_id <- paste0("q", o, "_staffing_five_levels")
    ppo <- .oe_capture(fit_ppo_extension(d, original, ppo_id))
    if (is.null(ppo$value)) {
      pf <- .failed_fit(ppo_id, paste("extension caught:", ppo$error), nrow(d), list(attempt = ppo))
      pf$G <- 8L; pf$type <- "ppo_ER"; pf$rows <- as.character(d$Observation)
    } else pf <- ppo$value
    factor_fit <- .oe_capture({
      df <- .oe_factor_data(d)
      if (any(table(df$B_f) == 0L) || any(table(df$A_f) == 0L)) stop("At least one original staffing category is absent.")
      f <- fit_ordinal(df, "Y", .oe_factor_predictors, model_id = fac_id)
      f$inference_valid <- FALSE; f$type <- "five_level_PO"; f$predictors <- .oe_factor_predictors
      # The original engine's ^(B_c|A_c) staffing-name check does not match
      # B_f/A_f. Apply its stricter 1e-4 tolerance to all eight indicators.
      factor_reasons <- oe_factor_alternate_reasons(f)
      f$diagnostics$factor_staffing_alternate_reasons <- factor_reasons
      if (length(factor_reasons)) {
        f$valid <- FALSE; f$status <- paste(c(f$status, factor_reasons), collapse = ";")
      }
      if (isTRUE(f$valid) && (length(f$coef) != 22L || f$logLik < original$logLik - 1e-6 * (1 + abs(original$logLik)))) {
        f$valid <- FALSE; f$status <- "factor dimension or nested likelihood failure"
      }
      if (isTRUE(f$valid)) f$status <- "accepted; point estimates only"
      f
    })
    if (is.null(factor_fit$value)) {
      ff <- .failed_fit(fac_id, paste("extension caught:", factor_fit$error), nrow(d), list(attempt = factor_fit))
      ff$G <- 8L; ff$type <- "five_level_PO"; ff$rows <- as.character(d$Observation)
    } else ff <- factor_fit$value
    for (f in list(pf, ff)) {
      fits_all[[f$model_id]] <- f; saveRDS(f, file.path(outdir, "fits", paste0(f$model_id, ".rds")))
      estimates[[f$model_id]] <- .oe_point_results(f, o)
      statuses[[f$model_id]] <- data.frame(model_id = f$model_id, outcome = o, n = f$n, sites = f$G,
        valid = f$valid, inference_valid = FALSE, status = f$status, log_likelihood = f$logLik)
    }
    support <- oe_staffing_support(d, o); supports[[as.character(o)]] <- support
    models <- setNames(list(original, pf, ff), c(original$model_id, pf$model_id, ff$model_id))
    predictions[[as.character(o)]] <- oe_probability_summary(models, d, o, support)
    saveRDS(list(support = support, predictions = predictions[[as.character(o)]]), file.path(outdir, paste0("q", o, "_probability_record.rds")))
    for (nm in c("status", "rows", "composition", "cells"))
      write.csv(support[[nm]], file.path(outdir, paste0("q", o, "_support_", nm, ".csv")), row.names = FALSE, na = "")
    for (nm in c("probabilities", "contrasts", "status"))
      write.csv(predictions[[as.character(o)]][[nm]], file.path(outdir, paste0("q", o, "_probability_", nm, ".csv")), row.names = FALSE, na = "")
  }
  write.csv(do.call(rbind, estimates), file.path(outdir, "ordinal_extension_point_estimates.csv"), row.names = FALSE, na = "")
  write.csv(do.call(rbind, statuses), file.path(outdir, "ordinal_extension_model_status.csv"), row.names = FALSE, na = "")
  write.csv(data.frame(path = source_files, md5 = unname(tools::md5sum(source_files))),
    file.path(outdir, "source_fit_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
  result <- list(fits = fits_all, estimates = do.call(rbind, estimates), status = do.call(rbind, statuses),
    support = supports, predictions = predictions, source_files = source_files,
    purpose = "Post-results specification 1.4 point-estimate sensitivity; no new inference")
  saveRDS(result, file.path(outdir, "ordinal_extensions.rds"))
  invisible(result)
}
