# Exact input decoding and prespecified analysis populations.
# No functions in this file fit models or execute on sourcing.
likert_labels <- c("Strongly disagree", "Disagree", "Neither agree nor disagree",
                   "Agree", "Strongly agree")
analysis_items <- paste0("Q", c(1:4, 6:9, 11:13, 15, 17, 21, 23:27))
domain_items <- list(L = paste0("Q", c(21, 23:27)), E = paste0("Q", 1:4),
                     Eminus = paste0("Q", 1:3), R = paste0("Q", 7:9),
                     M = paste0("Q", c(13, 15, 17)))

strict_factor <- function(x, levels, name) {
  bad <- unique(x[!is.na(x) & !x %in% levels])
  if (length(bad)) stop(name, ": undocumented labels: ", paste(bad, collapse = " | "))
  factor(x, levels = levels)
}

read_study_data <- function(path) {
  # Only the documented empty string is missing. Literal NA is an unknown label.
  raw <- read.csv(path, colClasses = "character", na.strings = "", check.names = FALSE,
                  strip.white = FALSE, stringsAsFactors = FALSE)
  expected <- c("Observation", "Site", "Language", paste0("Q", 1:28), paste0("D", 1:9))
  if (!identical(names(raw), expected)) stop("Input questionnaire columns or order differ from the expected schema.")
  if (!nrow(raw)) stop("Input has no rows.")
  if (anyNA(raw$Observation) || anyDuplicated(raw$Observation) ||
      any(!grepl("^[0-9]+$", raw$Observation))) stop("Observation identifiers must be unique nonmissing digit strings.")
  if (anyNA(raw$Site) || any(!raw$Site %in% as.character(1:8))) stop("Site must be observed and belong to 1:8.")
  if (!setequal(raw$Site, as.character(1:8))) stop("All eight expected sites must occur in the input.")
  dat <- raw
  for (q in paste0("Q", 1:28)) {
    strict_factor(raw[[q]], likert_labels, q)
    dat[[q]] <- match(raw[[q]], likert_labels)
  }
  tenure <- c("Less than 1 year", "1-5 years", "6-10 years", "11-15 years", "Greater than 15 years")
  dat$D1 <- strict_factor(raw$D1, tenure, "D1")
  dat$D5 <- strict_factor(raw$D5, tenure, "D5")
  dat$D2 <- strict_factor(raw$D2, c("false", "true"), "D2")
  dat$Language <- strict_factor(raw$Language, c("English", "Spanish"), "Language")
  dat$Site <- factor(raw$Site, levels = as.character(1:8))
  attr(dat, "raw_input") <- raw
  dat
}

required_predictors <- function(outcome) {
  stopifnot(outcome %in% c(4, 6))
  unique(c("Site", "Q11", "Q12", domain_items$L,
           if (outcome == 6) domain_items$E else domain_items$Eminus,
           domain_items$R, domain_items$M))
}

complete_case_rows <- function(data, outcome) {
  complete.cases(data[, unique(c(paste0("Q", outcome), required_predictors(outcome))), drop = FALSE])
}

derive_scores <- function(data, outcome = 6, management = "stored", outcome_coding = "four") {
  stopifnot(outcome %in% c(4, 6), management %in% c("stored", "drop_Q13", "reverse_Q13", "reverse_Q15", "reverse_Q17", "reverse_all"))
  d <- data
  d$Site <- factor(as.character(d$Site), levels = as.character(1:8))
  items <- domain_items
  if (management == "drop_Q13") items$M <- c("Q15", "Q17")
  management_data <- d[, items$M, drop = FALSE]
  reverse <- switch(management, reverse_Q13 = "Q13", reverse_Q15 = "Q15", reverse_Q17 = "Q17",
                    reverse_all = items$M, character())
  for (q in reverse) management_data[[q]] <- 6 - management_data[[q]]
  for (nm in c("L", "E", "Eminus", "R")) d[[nm]] <- rowMeans(d[, items[[nm]], drop = FALSE], na.rm = FALSE)
  d$M <- rowMeans(management_data, na.rm = FALSE)
  d$B <- d$Q11
  d$A <- d$Q12
  for (nm in c("L", "E", "Eminus", "R", "M", "B", "A")) d[[paste0(nm, "_c")]] <- d[[nm]] - 3
  if (outcome == 4) d$E_c <- d$Eminus_c
  y <- d[[paste0("Q", outcome)]]
  if (outcome_coding == "four") {
    d$Y <- ordered(ifelse(is.na(y), NA_integer_, pmax(y - 1L, 1L)), levels = 1:4)
  } else if (outcome_coding == "five" && outcome == 6) {
    d$Y <- ordered(y, levels = 1:5)
  } else if (outcome_coding == "three" && outcome == 4) {
    d$Y <- ordered(ifelse(is.na(y), NA_integer_, pmax(y - 2L, 1L)), levels = 1:3)
  } else stop("Unspecified outcome coding.")
  d
}

staffing_factor <- function(x) {
  factor(ifelse(is.na(x), NA_character_, ifelse(x <= 2, "disagreement", ifelse(x == 3, "neutral", "agreement"))),
         levels = c("neutral", "disagreement", "agreement"))
}

outcome_distribution <- function(x, levels, population, outcome, coding) {
  n <- tabulate(match(x, levels), nbins = length(levels))
  denom <- sum(!is.na(x))
  data.frame(population = population, outcome = outcome, coding = coding,
             category = levels, n = n, denominator_observed = denom,
             percent = if (denom) 100 * n / denom else NA_real_, missing = sum(is.na(x)))
}

write_input_descriptives <- function(data, outdir, source_path) {
  raw <- attr(data, "raw_input")
  if (is.null(raw)) stop("Input manifest requires the original decoded input attribute.")
  write.csv(data.frame(input = normalizePath(source_path, winslash = "/"), rows = nrow(data), columns = ncol(raw),
                       missing_token = "empty string only", management_direction = "as stored; transformation provenance unconfirmed"),
            file.path(outdir, "input_manifest.csv"), row.names = FALSE)
  freq <- do.call(rbind, lapply(names(raw), function(nm) {
    tt <- table(raw[[nm]], useNA = "always")
    data.frame(variable = nm, stored_label = names(tt), n = as.integer(tt))
  }))
  write.csv(freq, file.path(outdir, "input_category_frequencies.csv"), row.names = FALSE, na = "")
  mapping <- data.frame(stored_label = likert_labels, numeric_code = 1:5)
  write.csv(mapping, file.path(outdir, "likert_mapping.csv"), row.names = FALSE)
  variables <- unique(c("Site", analysis_items, "D1", "D2", "D5", "Language"))
  ds <- derive_scores(data)
  variables <- c(variables, "B", "A", "L", "E", "Eminus", "R", "M")
  missingness <- do.call(rbind, lapply(c("all", as.character(1:8)), function(site) {
    idx <- if (site == "all") rep(TRUE, nrow(ds)) else as.character(ds$Site) == site
    do.call(rbind, lapply(variables, function(v) {
      data.frame(site = site, variable = v, n = sum(idx), missing = sum(is.na(ds[[v]][idx])),
                 percent_missing = 100 * mean(is.na(ds[[v]][idx])))
    }))
  }))
  write.csv(missingness, file.path(outdir, "missingness_overall_by_site.csv"), row.names = FALSE)
  distributions <- list(); flow <- list(); exclusion <- list()
  for (o in c(6, 4)) {
    d <- derive_scores(data, o)
    pred_complete <- complete.cases(data[, required_predictors(o), drop = FALSE])
    cc <- complete_case_rows(data, o)
    outcome_seen <- !is.na(data[[paste0("Q", o)]])
    flow[[as.character(o)]] <- data.frame(outcome = o, population = c("input", "outcome_observed", "predictors_complete", "outcome_and_predictors_complete", "worker_covariates_complete"),
        n = c(nrow(data), sum(outcome_seen), sum(pred_complete), sum(cc), sum(cc & complete.cases(data[, c("D1", "D2", "D5", "Language")]))))
    exclusion[[as.character(o)]] <- data.frame(Observation = data$Observation, outcome = o, included = cc,
      missing_required = apply(is.na(data[, unique(c(paste0("Q", o), required_predictors(o))), drop = FALSE]), 1,
                              function(z) paste(names(z)[z], collapse = ";")))
    groups <- list(all = rep(TRUE, nrow(data)), predictors_complete = pred_complete,
                   predictors_incomplete = !pred_complete, complete_case = cc)
    for (g in names(groups)) {
      idx <- groups[[g]]
      distributions[[paste(o, g, 5)]] <- outcome_distribution(data[[paste0("Q", o)]][idx], 1:5, g, o, "five")
      distributions[[paste(o, g, 4)]] <- outcome_distribution(as.integer(d$Y[idx]), 1:4, g, o, "four")
    }
  }
  write.csv(do.call(rbind, flow), file.path(outdir, "participant_flow.csv"), row.names = FALSE)
  write.csv(do.call(rbind, exclusion), file.path(outdir, "participant_exclusions.csv"), row.names = FALSE)
  write.csv(do.call(rbind, distributions), file.path(outdir, "outcome_distributions.csv"), row.names = FALSE)
  ie <- complete.cases(ds[, c("Site", domain_items$L, "Q12", "Q6"), drop = FALSE])
  e <- ds[ie, , drop = FALSE]
  e$leadership_high <- e$L >= 4
  e$absence_category <- ifelse(e$A <= 2, "inadequate", ifelse(e$A == 3, "neutral", "adequate"))
  cells <- expand.grid(leadership_high = c(FALSE, TRUE), absence_category = c("inadequate", "neutral", "adequate"), stringsAsFactors = FALSE)
  rq4 <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    idx <- e$leadership_high == cells$leadership_high[i] & e$absence_category == cells$absence_category[i]
    n <- sum(idx)
    data.frame(leadership_high = cells$leadership_high[i], absence_category = cells$absence_category[i],
      n = n, denominator = nrow(e), percent = if (nrow(e)) 100 * n / nrow(e) else NA_real_,
      sparse = n < 10, outcome_category = 1:4, outcome_n = tabulate(as.integer(e$Y[idx]), nbins = 4),
      outcome_percent = if (n) 100 * tabulate(as.integer(e$Y[idx]), nbins = 4) / n else rep(NA_real_, 4))
  }))
  write.csv(rq4, file.path(outdir, "rq4_six_cells.csv"), row.names = FALSE)
  joint <- as.data.frame(table(Q11 = factor(data$Q11, levels = 1:5), Q12 = factor(data$Q12, levels = 1:5), useNA = "always"))
  write.csv(joint, file.path(outdir, "staffing_joint_frequencies_all.csv"), row.names = FALSE)
  staff_observed <- complete.cases(data[, c("Q11", "Q12"), drop = FALSE])
  b <- data$Q11[staff_observed]; a <- data$Q12[staff_observed]
  variable <- length(b) > 1L && length(unique(b)) > 1L && length(unique(a)) > 1L
  write.csv(data.frame(population = "both_staffing_items_observed", n = sum(staff_observed),
    pearson = if (variable) cor(b, a, method = "pearson") else NA_real_,
    spearman = if (variable) cor(b, a, method = "spearman") else NA_real_),
    file.path(outdir, "staffing_correlations_all.csv"), row.names = FALSE)
  invisible(TRUE)
}
