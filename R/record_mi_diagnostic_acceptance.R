# Record a completed human/agent diagnostic inspection after structural checks.
# This script does not diagnose convergence automatically or fit outcome models.
# Rscript record_mi_diagnostic_acceptance.R RUN_DIRECTORY REVIEW_MARKDOWN REVIEWER_NAME
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L, nzchar(trimws(args[3])))
run_dir <- args[1]; review_file <- args[2]
state <- readRDS(file.path(run_dir, "mi_checkpoint.rds"))
stopifnot(state$status == "diagnostics_pending", file.exists(review_file),
          identical(state$seeds, c(20261111L, 20261112L)))
items <- paste0("Q", c(1:4,6:9,11:13,15,17,21,23:27))
checks <- list()
for (seed in state$seeds) {
  chains <- state$runs[[as.character(seed)]]
  stopifnot(length(chains) == state$target_m,
            identical(vapply(chains, function(x) x$chain, integer(1)), seq_len(state$target_m)))
  for (chain in chains) {
    stopifnot(chain$status == "complete", chain$iteration == state$target_iterations,
              identical(chain$data$Observation, state$original$Observation),
              identical(chain$data$Site, state$original$Site),
              length(chain$trace) == state$target_iterations)
    for (item in items) {
      observed <- !is.na(state$original[[item]])
      stopifnot(identical(chain$data[[item]][observed], state$original[[item]][observed]),
                all(chain$data[[item]] %in% 1:5))
    }
    checks[[length(checks)+1L]] <- data.frame(seed=seed, chain=chain$chain,
      iterations=chain$iteration, original_entries_preserved=TRUE, support_valid=TRUE)
  }
}
stem <- paste0("review_m",state$target_m,"_iter",state$target_iterations)
write.csv(do.call(rbind,checks),file.path(run_dir,paste0(stem,"_structural_checks.csv")),row.names=FALSE)
review <- list(reviewer=args[3],notes=paste(readLines(review_file,warn=FALSE),collapse="\n"),decision="accept")
saveRDS(review,file.path(run_dir,paste0(stem,".rds")))
cat("Accepted structural verification for",length(checks),"completed chains; visual review retained in",review_file,"\n")
