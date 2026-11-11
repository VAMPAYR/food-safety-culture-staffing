# Synthetic input and population verification. Sourcing performs no execution.
run_data_prep_tests <- function(prep_path, output_dir = NULL) {
  source(prep_path, local = environment())
  results <- list()
  check <- function(name, pass, detail = "") {
    results[[length(results) + 1L]] <<- data.frame(test = name, passed = isTRUE(pass),
      detail = as.character(detail), stringsAsFactors = FALSE)
  }
  rejects <- function(expr, pattern) {
    problem <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
    !is.null(problem) && grepl(pattern, problem)
  }
  n <- 16L
  base <- data.frame(Observation = paste0("100", seq_len(n)),
                     Site = factor(rep(1:8, each = 2L), levels = 1:8),
                     Language = factor(rep(c("English", "Spanish"), 8L), levels = c("English", "Spanish")))
  for (q in 1:28) base[[paste0("Q", q)]] <- ((seq_len(n) + q - 2L) %% 5L) + 1L
  tenure <- c("Less than 1 year", "1-5 years", "6-10 years", "11-15 years", "Greater than 15 years")
  for (j in 1:9) base[[paste0("D", j)]] <- rep("unused", n)
  base$D1 <- factor(rep(tenure, length.out = n), levels = tenure)
  base$D5 <- factor(rep(rev(tenure), length.out = n), levels = tenure)
  base$D2 <- factor(rep(c("false", "true"), 8L), levels = c("false", "true"))
  check("synthetic_schema_28_questionnaire_19_analysis_items",
    length(analysis_items) == 19L && all(paste0("Q", 1:28) %in% names(base)))

  d4 <- derive_scores(base, 4)
  d6 <- derive_scores(base, 6)
  direct_minus <- (base$Q1 + base$Q2 + base$Q3) / 3
  direct_full <- (base$Q1 + base$Q2 + base$Q3 + base$Q4) / 4
  check("Q4_employee_score_omits_outcome", max(abs(d4$E_c - (direct_minus - 3))) < 1e-12 &&
    !"Q4" %in% required_predictors(4) && "Q4" %in% required_predictors(6))
  check("Q6_employee_score_keeps_Q4", max(abs(d6$E_c - (direct_full - 3))) < 1e-12)
  changed <- base
  changed$Q4[1] <- if (base$Q4[1] == 5L) 1L else 5L
  changed4 <- derive_scores(changed, 4)
  check("changing_Q4_changes_outcome_but_not_Q4_employee_predictor",
    changed4$Y[1] != d4$Y[1] && identical(changed4$E_c, d4$E_c))
  check("scores_center_on_three_without_standardization", all(d6$B_c == base$Q11 - 3) &&
    all(d6$A_c == base$Q12 - 3) && all(d6$L_c == d6$L - 3) && all(d6$M_c == d6$M - 3))

  missing <- base
  missing$Q4[2] <- NA_integer_
  missing$Q6[3] <- NA_integer_
  missing$Q21[4] <- NA_integer_
  missing$Q11[5] <- NA_integer_
  missing$Q15[6] <- NA_integer_
  missing$Q13[7] <- NA_integer_
  missing$D1[8] <- NA
  missing$D5[9] <- NA
  missing$D2[10] <- NA
  missing$Language[11] <- NA
  cc6 <- complete_case_rows(missing, 6)
  cc4 <- complete_case_rows(missing, 4)
  check("Q6_complete_cases_require_outcome_staffing_and_all_domain_items",
    identical(which(!cc6), 2:7))
  check("Q4_population_does_not_require_Q6",
    identical(which(!cc4), c(2L, 4L, 5L, 6L, 7L)) && cc4[3])
  p6 <- complete.cases(missing[, required_predictors(6), drop = FALSE])
  p4 <- complete.cases(missing[, required_predictors(4), drop = FALSE])
  check("predictor_completeness_separate_from_outcome_observation",
    p6[3] && !cc6[3] && p4[2] && !cc4[2] && p4[3])
  masks <- list(Q6 = !is.na(missing$Q6), Q4 = !is.na(missing$Q4))
  completed <- missing
  for (q in analysis_items) completed[[q]][is.na(completed[[q]])] <- 3L
  check("original_outcome_masks_survive_auxiliary_completion",
    identical(which(!masks$Q6), 3L) && identical(which(!masks$Q4), 2L) &&
      !anyNA(completed$Q6) && !anyNA(completed$Q4) &&
      !missing$Observation[3] %in% completed$Observation[masks$Q6] &&
      !missing$Observation[2] %in% completed$Observation[masks$Q4])
  dm <- derive_scores(missing, 6)
  check("domain_means_require_every_constituent_item",
    is.na(dm$E[2]) && is.na(dm$L[4]) && is.na(dm$M[6]) && is.na(dm$M[7]) &&
      !is.na(dm$Eminus[2]) && is.na(dm$B_c[5]))

  for (outcome in c(6L, 4L)) {
    idx <- complete_case_rows(missing, outcome)
    fixed <- missing[idx, , drop = FALSE]
    main <- derive_scores(fixed, outcome)
    variants <- c("stored", "drop_Q13", "reverse_Q13", "reverse_Q15", "reverse_Q17")
    for (variant in variants) {
      alt <- derive_scores(fixed, outcome, management = variant)
      check(paste0("Q", outcome, "_", variant, "_retains_original_complete_case_rows"),
        identical(alt$Observation, main$Observation) && nrow(alt) == sum(idx) &&
          !anyNA(alt[, c("Y", "L_c", "E_c", "R_c", "M_c", "B_c", "A_c")]))
    }
    alt <- derive_scores(fixed, outcome, outcome_coding = if (outcome == 6L) "five" else "three")
    check(paste0("Q", outcome, "_alternative_outcome_retains_original_rows"),
      identical(alt$Observation, main$Observation))
    check(paste0("Q", outcome, "_site_has_fixed_eight_level_schema"),
      identical(levels(main$Site), as.character(1:8)))
    subset <- fixed[complete.cases(fixed[, c("D1", "D5", "D2", "Language")]), , drop = FALSE]
    expected_worker_ids <- missing$Observation[idx & !seq_len(n) %in% 8:11]
    check(paste0("Q", outcome, "_worker_subset_rows_identical_for_base_and_augmented_designs"),
      identical(derive_scores(subset, outcome)$Observation, expected_worker_ids) &&
        identical(subset$Observation, expected_worker_ids) && nrow(subset) == sum(idx) - 4L)
  }
  dropped <- derive_scores(missing, 6, management = "drop_Q13")
  check("dropping_Q13_must_not_expand_original_cohort", !is.na(dropped$M[7]) && !cc6[7])
  check("management_drop_Q13_score", all(derive_scores(base, 6, "drop_Q13")$M == (base$Q15 + base$Q17) / 2))
  for (q in c("Q13", "Q15", "Q17")) {
    expected <- (base$Q13 + base$Q15 + base$Q17 + 6 - 2 * base[[q]]) / 3
    check(paste0("selective_", q, "_reversal_score"),
      max(abs(derive_scores(base, 6, paste0("reverse_", q))$M - expected)) < 1e-12)
  }
  for (o in c(6L, 4L)) {
    codes <- base
    codes[[paste0("Q", o)]][1:5] <- 1:5
    four <- derive_scores(codes, o)
    check(paste0("Q", o, "_four_category_mapping_and_order"),
      is.ordered(four$Y) && identical(levels(four$Y), as.character(1:4)) &&
        identical(as.integer(four$Y[1:5]), c(1L, 1L, 2L, 3L, 4L)))
    alternative <- derive_scores(codes, o, outcome_coding = if (o == 6L) "five" else "three")
    expected <- if (o == 6L) 1:5 else c(1L, 1L, 1L, 2L, 3L)
    check(paste0("Q", o, "_alternative_category_mapping"), identical(as.integer(alternative$Y[1:5]), expected))
  }
  staff <- staffing_factor(c(1:5, NA_integer_))
  check("staffing_factor_neutral_reference_and_three_levels",
    identical(levels(staff), c("neutral", "disagreement", "agreement")) &&
      identical(as.character(staff), c("disagreement", "disagreement", "neutral", "agreement", "agreement", NA_character_)))

  # Exercise the actual strict CSV reader using synthetic labels only.
  fixture_dir <- if (is.null(output_dir)) tempfile("data_prep_fixtures_") else file.path(output_dir, "synthetic_input_fixtures")
  dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)
  raw <- base
  for (q in paste0("Q", 1:28)) raw[[q]] <- likert_labels[base[[q]]]
  fixture <- file.path(fixture_dir, "valid_synthetic.csv")
  write.csv(raw, fixture, row.names = FALSE, na = "")
  decoded <- read_study_data(fixture)
  check("exact_stored_labels_decode_without_recoding", all(vapply(paste0("Q", 1:28),
    function(q) identical(decoded[[q]], base[[q]]), logical(1))))
  check("worker_covariates_retain_full_categories_and_references",
    identical(levels(decoded$D1), tenure) && identical(levels(decoded$D5), tenure) &&
      identical(levels(decoded$D2), c("false", "true")) &&
      identical(levels(decoded$Language), c("English", "Spanish")))
  for (token in c("NA", " Agree", "unknown")) {
    invalid <- raw; invalid$Q6[1] <- token
    path <- file.path(fixture_dir, paste0("invalid_label_", match(token, c("NA", " Agree", "unknown")), ".csv"))
    write.csv(invalid, path, row.names = FALSE, na = "")
    check(paste0("unknown_label_rejected_", token), rejects(read_study_data(path), "undocumented labels"))
  }
  blank <- raw; blank$Q6[1] <- ""
  path <- file.path(fixture_dir, "blank_synthetic.csv")
  write.csv(blank, path, row.names = FALSE, na = "")
  check("documented_blank_only_is_missing", is.na(read_study_data(path)$Q6[1]))
  duplicate <- raw; duplicate$Observation[2] <- duplicate$Observation[1]
  path <- file.path(fixture_dir, "duplicate_synthetic.csv")
  write.csv(duplicate, path, row.names = FALSE, na = "")
  check("duplicate_identifier_rejected", rejects(read_study_data(path), "identifiers"))
  invalid_site <- raw; invalid_site$Site <- as.character(invalid_site$Site); invalid_site$Site[1] <- "9"
  path <- file.path(fixture_dir, "invalid_site_synthetic.csv")
  write.csv(invalid_site, path, row.names = FALSE, na = "")
  check("unknown_site_rejected", rejects(read_study_data(path), "Site"))

  out <- do.call(rbind, results)
  if (!is.null(output_dir)) write.csv(out, file.path(output_dir, "data_prep_tests.csv"), row.names = FALSE)
  if (any(!out$passed)) stop("Data preparation verification failed: ", paste(out$test[!out$passed], collapse = ", "))
  out
}
