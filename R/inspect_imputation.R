# Diagnostic presentation only. Reads saved traces; does not generate data or fit models.
# Rscript inspect_imputation.R --run-dir OUTPUT_DIRECTORY
inspect_imputation <- function(run_dir) {
  trace <- read.csv(file.path(run_dir, "mi_trace.csv"))
  chains <- read.csv(file.path(run_dir, "mi_chains.csv"))
  if (!nrow(trace) || any(chains$status != "complete")) stop("Complete chain diagnostics required.")
  out <- file.path(run_dir, paste0("diagnostic_review_m", length(unique(chains$chain)), "_iter", max(chains$iteration)))
  dir.create(out, showWarnings = FALSE)
  items <- unique(trace$item)[order(as.integer(sub("Q", "", unique(trace$item))))]
  summaries <- list()
  for (seed in unique(trace$seed)) {
    for (moment in c("mean", "variance")) {
      chunks <- split(items, ceiling(seq_along(items) / 9))
      for (page in seq_along(chunks)) {
        png(file.path(out, paste0("trace_", seed, "_", moment, "_", page, ".png")), width = 1350, height = 1200, res = 140)
        par(mfrow = c(3, 3), mar = c(3, 3, 2.4, 0.5), oma = c(0, 0, 2, 0))
        for (item in chunks[[page]]) {
          z <- trace[trace$seed == seed & trace$item == item, ]
          finite <- is.finite(z[[moment]])
          if (!any(finite)) {
            plot.new(); title(paste(item, "n missing", z$n_imputed[1])); text(0.5, 0.5, "Variance undefined for one entry")
            next
          }
          yr <- range(z[[moment]][finite]); if (diff(yr) == 0) yr <- yr + c(-0.5, 0.5)
          plot(range(z$iteration), yr, type = "n", xlab = "Iteration", ylab = moment,
               main = paste(item, "n missing", z$n_imputed[1]))
          for (chain in unique(z$chain)) {
            w <- z[z$chain == chain, ]; lines(w$iteration, w[[moment]], col = adjustcolor("gray40", alpha.f = 0.18))
          }
          average <- aggregate(z[[moment]], list(iteration = z$iteration), mean, na.rm = TRUE)
          lines(average$iteration, average$x, col = "#9b2525", lwd = 2)
        }
        mtext(paste("Root seed", seed, "Imputed entry", moment, "Red line is the mean across chains"), outer = TRUE, cex = 1)
        dev.off()
      }
    }
    for (item in items) {
      z <- trace[trace$seed == seed & trace$item == item, ]
      last <- max(z$iteration); span <- min(10L, floor(last/2))
      early <- z[z$iteration <= span, ]; late <- z[z$iteration > last-span, ]
      summaries[[length(summaries)+1L]] <- data.frame(seed = seed, item = item,
        n_imputed = z$n_imputed[1], iterations = last,
        early_mean = mean(early$mean), late_mean = mean(late$mean),
        difference = mean(late$mean)-mean(early$mean),
        early_variance = mean(early$variance, na.rm = TRUE), late_variance = mean(late$variance, na.rm = TRUE))
    }
  }
  write.csv(do.call(rbind, summaries), file.path(out, "trace_early_late_summary.csv"), row.names = FALSE)
  distributions <- read.csv(file.path(run_dir, "mi_distributions.csv"))
  averaged <- aggregate(distributions[, c("n", "denominator")], distributions[, c("seed", "item", "site", "kind", "value")], mean)
  averaged$percent <- ifelse(averaged$denominator > 0, 100 * averaged$n / averaged$denominator, NA_real_)
  write.csv(averaged, file.path(out, "mean_distributions_by_site.csv"), row.names = FALSE)
  donors <- read.csv(file.path(run_dir, "mi_donor_reuse.csv"))
  reuse <- do.call(rbind, lapply(split(donors, interaction(donors$seed, donors$item, drop = TRUE)), function(z) {
    data.frame(seed = z$seed[1], item = z$item[1],
      distinct_donors_terminal = length(unique(z$donor_row[z$uses_terminal_iteration > 0])),
      maximum_uses_in_one_terminal_chain = max(z$uses_terminal_iteration),
      maximum_uses_in_one_chain_all_iterations = max(z$uses_all_iterations))
  }))
  write.csv(reuse, file.path(out, "donor_reuse_summary.csv"), row.names = FALSE)
  correlations <- read.csv(file.path(run_dir, "mi_correlations.csv"))
  summaries <- do.call(rbind, lapply(split(correlations, interaction(correlations$seed, correlations$scope, correlations$item_a, correlations$item_b, drop = TRUE)), function(z) {
    values <- z$correlation[is.finite(z$correlation)]
    data.frame(seed = z$seed[1], scope = z$scope[1], item_a = z$item_a[1], item_b = z$item_b[1],
      n = z$n[1], valid_chains = length(values),
      mean_correlation = if (length(values)) mean(values) else NA_real_,
      minimum = if (length(values)) min(values) else NA_real_, maximum = if (length(values)) max(values) else NA_real_)
  }))
  write.csv(summaries, file.path(out, "correlation_summary.csv"), row.names = FALSE)
  writeLines(c("These summaries support visual diagnostic inspection; they are not hypothesis tests.",
    "No numeric drift cutoff is imposed. The recorded reviewer must consider persistent drift, mixing, support, and donor concentration.",
    "Few missing entries can produce volatile means or undefined variances and correlations."), file.path(out, "README.txt"))
  message("Diagnostic review files: ", normalizePath(out, winslash = "/"))
  invisible(out)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L || args[1] != "--run-dir") stop("Use --run-dir OUTPUT_DIRECTORY")
  inspect_imputation(args[2])
}
