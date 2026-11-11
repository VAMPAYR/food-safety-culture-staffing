# Native Windows CLI for specification 1.4 priority 3. Definition-only on source.
# Example (after independent mathematical and implementation review):
# Rscript R/run_imputation_extensions.R --stage generate --output outputs/public_priority3_20261111 --m 50 --iterations 50 --workers 4
# Inspect diagnostics, then record-review --decision accept --reviewer NAME --note TEXT.
# Fit is a separate invocation with --stage fit. No generation call fits outcomes.
mie_sha256 <- function(path) {
  executable <- file.path(Sys.getenv("SystemRoot","C:/Windows"),"System32","certutil.exe")
  if (!file.exists(executable)) stop("Windows certutil is required for SHA256 provenance.")
  output <- system2(executable,c("-hashfile",shQuote(normalizePath(path,winslash="\\",mustWork=TRUE)),"SHA256"),stdout=TRUE,stderr=TRUE)
  hash <- trimws(output[grepl("^[[:xdigit:]]{64}$",trimws(output))])
  if(length(hash)!=1L) stop("SHA256 failed: ",path)
  toupper(hash)
}
mie_cli_options <- function(args) {
  options <- list(stage="generate",output=NULL,data=NULL,m=50L,iterations=50L,workers=4L,scenario="all",decision=NULL,reviewer=NULL,note=NULL)
  if(length(args)%%2L) stop("Use --name value pairs.")
  if(length(args)) for(i in seq(1L,length(args),2L)) {
    name <- sub("^--","",args[i]); if(!startsWith(args[i],"--") || !name %in% names(options)) stop("Unknown option: ",args[i])
    options[[name]] <- args[i+1L]
  }
  if(!options$stage %in% c("generate","diagnostics","record-review","fit","summary")) stop("Unknown stage.")
  if(is.null(options$output)) stop("--output is required.")
  for(name in c("m","iterations","workers")) options[[name]] <- as.integer(options[[name]])
  if(!options$m %in% c(50L,100L,200L) || !options$iterations %in% c(50L,100L,200L) || !options$workers %in% 1:4) stop("Invalid chain/iteration/worker count.")
  options
}
mie_source_paths <- function(script_dir) {
  # Public clones omit the internal mathematical specification. Every executable
  # source remains mandatory; local runs additionally bind the specification.
  required <- file.path(script_dir,c("data_prep.R","model_engine.R","imputation_extensions.R","run_imputation_extensions.R","test_imputation_extensions.R"))
  if(any(!file.exists(required))) stop("A required reviewed R source is absent.")
  specification <- file.path(dirname(script_dir),"specification","mathematical_specification","main.tex")
  c(required,if(file.exists(specification))specification else character())
}
run_imputation_scenario <- function(args=commandArgs(trailingOnly=TRUE),script_dir=NULL) {
  if(is.null(script_dir)) {
    script <- grep("^--file=",commandArgs(FALSE),value=TRUE)
    if(length(script)!=1L) stop("Supply script_dir when called interactively.")
    script_dir <- dirname(normalizePath(sub("^--file=","",script)))
  }
  environment <- environment(run_imputation_scenario)
  for(file in c("data_prep.R","model_engine.R","imputation_extensions.R")) source(file.path(script_dir,file),local=environment)
  options <- mie_cli_options(args); outdir <- options$output
  if(!options$scenario %in% mie_scenarios()$scenario) stop("Select one recognized scenario.")
  directory_exists <- dir.exists(outdir)
  if(!directory_exists && options$stage!="generate") stop("An existing generation directory is required.")
  if(!directory_exists) dir.create(outdir,recursive=TRUE)
  outdir <- normalizePath(outdir,winslash="/",mustWork=TRUE)
  manifest_path <- file.path(outdir,"manifest.rds")
  source_files <- mie_source_paths(script_dir)
  source_hashes <- data.frame(file=normalizePath(source_files,winslash="/"),sha256=vapply(source_files,mie_sha256,character(1L)),row.names=NULL)
  if(!file.exists(manifest_path)) {
    if(directory_exists && length(list.files(outdir,all.files=TRUE,no..=TRUE))) stop("New generation output must be empty.")
    if(options$m!=50L || options$iterations!=50L) stop("Initial generation requires 50 chains and 50 iterations.")
    input <- if(is.null(options$data))file.path(dirname(dirname(script_dir)),"dataset","data.csv") else options$data
    input_hash <- mie_sha256(input)
    if(input_hash!="2DFCE56408F4411C01426556C250DBD932790CFD24293B278E63F32A8323EDED") stop("Raw dataset digest does not match reviewed input.")
    original <- mie_encode(read_study_data(input)); mie_validate(original)
    manifest <- list(version="1.4",input=normalizePath(input,winslash="/"),input_sha256=input_hash,source_hashes=source_hashes,
      m=50L,iterations=50L,seeds=mie_seeds(),scenario=options$scenario,created=as.character(Sys.time()),post_results_amendment=TRUE)
    mie_save(original,file.path(outdir,"original_encoded.rds")); manifest$original_md5 <- unname(tools::md5sum(file.path(outdir,"original_encoded.rds"))); mie_save(manifest,manifest_path)
    write.csv(source_hashes,file.path(outdir,"source_sha256.csv"),row.names=FALSE)
    writeLines(capture.output(sessionInfo()),file.path(outdir,"sessionInfo.txt"))
  } else {
    manifest <- readRDS(manifest_path); original <- readRDS(file.path(outdir,"original_encoded.rds")); mie_validate(original)
    if(manifest$scenario!=options$scenario || !identical(unname(tools::md5sum(file.path(outdir,"original_encoded.rds"))),manifest$original_md5)) stop("Original data or scenario differs from manifest.")
    if(!identical(source_hashes,manifest$source_hashes)) stop("Code/specification changed since generation; use a new reviewed output directory.")
    if(!is.null(options$data) && mie_sha256(options$data)!=manifest$input_sha256) stop("Resume input digest differs.")
  }
  review_path <- function(m,iterations) file.path(outdir,sprintf("review_m%d_i%d.rds",m,iterations))
  if(options$stage=="generate") {
    changed_m <- options$m!=manifest$m; changed_iterations <- options$iterations!=manifest$iterations
    if(changed_m && changed_iterations) stop("Extend chains and iterations in separate reviewed steps.")
    if(changed_iterations) {
      next_iterations <- if(manifest$iterations==50L)100L else if(manifest$iterations==100L)200L else NA_integer_
      review <- if(file.exists(review_path(manifest$m,manifest$iterations)))readRDS(review_path(manifest$m,manifest$iterations)) else NULL
      if(is.null(review)||review$decision!="extend_iterations"||is.na(next_iterations)||options$iterations!=next_iterations) stop("Iteration extension needs the recorded next-step diagnostic decision.")
      previous_index <- readRDS(file.path(outdir,sprintf("diagnostic_index_m%d_i%d.rds",manifest$m,manifest$iterations)))
      mie_require_review(manifest_path,previous_index,review_path(manifest$m,manifest$iterations),require_accept=FALSE)
    }
    if(changed_m) {
      next_m <- if(manifest$m==50L)100L else if(manifest$m==100L)200L else NA_integer_
      previous_summary <- file.path(outdir,sprintf("summaries_m%d_i%d",manifest$m,manifest$iterations),"summary.rds")
      if(!file.exists(previous_summary)||!isTRUE(readRDS(previous_summary)$extension_required)||is.na(next_m)||options$m!=next_m) stop("Chain extension needs the specified Monte Carlo criterion and next count.")
    }
    manifest$m <- options$m; manifest$iterations <- options$iterations; mie_save(manifest,manifest_path)
    streams <- setNames(lapply(mie_seeds(),function(seed)mie_streams(seed,manifest$m)),as.character(mie_seeds()))
    tasks <- mie_tasks(outdir,manifest$m,manifest$iterations,manifest$scenario)
    results <- mie_parallel(tasks,mie_generation_worker,script_dir,original,streams,options$workers)
    write.csv(mie_bind(results),file.path(outdir,sprintf("generation_registry_m%d_i%d.csv",manifest$m,manifest$iterations)),row.names=FALSE)
    diagnostics <- mie_write_diagnostics(outdir,original,manifest$m,manifest$iterations,manifest$scenario)
    mie_save(diagnostics,file.path(outdir,sprintf("diagnostic_index_m%d_i%d.rds",manifest$m,manifest$iterations)))
    message("Generation stopped at diagnostic review gate: ",diagnostics$directory)
    return(invisible(diagnostics))
  }
  if(options$m!=manifest$m || options$iterations!=manifest$iterations) stop("Requested counts differ from current manifest.")
  index_path <- file.path(outdir,sprintf("diagnostic_index_m%d_i%d.rds",manifest$m,manifest$iterations))
  if(options$stage=="diagnostics") {
    diagnostics <- mie_write_diagnostics(outdir,original,manifest$m,manifest$iterations,manifest$scenario); mie_save(diagnostics,index_path)
    return(invisible(diagnostics))
  }
  if(!file.exists(index_path)) stop("Generate current diagnostics first.")
  diagnostics <- readRDS(index_path)
  if(options$stage=="record-review") {
    if(is.null(options$decision)||!options$decision %in% c("accept","extend_iterations","unavailable") || is.null(options$reviewer)||!nzchar(options$reviewer)||is.null(options$note)||!nzchar(options$note)) stop("Review requires --decision, --reviewer, and --note.")
    mie_validate_diagnostics(diagnostics)
    if(options$decision=="accept" && !isTRUE(diagnostics$complete)) stop("Incomplete or failed chains cannot be accepted.")
    if(options$decision=="extend_iterations" && manifest$iterations==200L) stop("At iteration cap; record unavailable.")
    if(options$decision=="unavailable" && isTRUE(diagnostics$complete) && manifest$iterations<200L) stop("Unresolved mixing/drift must extend to 100 then 200 iterations before unavailable.")
    path <- review_path(manifest$m,manifest$iterations)
    if(file.exists(path)) stop("Review already recorded; do not overwrite its decision.")
    review <- list(decision=options$decision,reviewer=options$reviewer,note=options$note,time=as.character(Sys.time()),
      manifest_md5=unname(tools::md5sum(manifest_path)),diagnostic_fingerprint=diagnostics$fingerprint)
    mie_save(review,path); return(invisible(review))
  }
  if(options$stage=="fit") {
    path <- review_path(manifest$m,manifest$iterations)
    if(!file.exists(path)) stop("A recorded diagnostic review is required before outcome fitting.")
    review <- mie_require_review(manifest_path,diagnostics,path)
    tasks <- mie_fit_tasks(outdir,manifest$m,manifest$iterations,manifest$scenario)
    registry <- mie_parallel(tasks,mie_fit_worker,script_dir,original,workers=options$workers)
  } else {
    tasks <- mie_fit_tasks(outdir,manifest$m,manifest$iterations,manifest$scenario)
    registry <- list()
    for(i in seq_len(nrow(tasks))) {
      fit <- mie_load_saved_fit(as.list(tasks[i,,drop=FALSE]))
      if(is.null(fit)) next
      registry[[length(registry)+1L]] <- fit
    }
  }
  summaries <- mie_write_fit_summary(registry,outdir,manifest$m,manifest$iterations,manifest$scenario)
  model_files <- c(tasks$fit_path,paste0(tasks$fit_path,".failure.rds")); model_files <- model_files[file.exists(model_files)]
  if(length(model_files)) write.csv(data.frame(file=model_files,md5=unname(tools::md5sum(model_files))),
    file.path(outdir,sprintf("model_md5_m%d_i%d.csv",manifest$m,manifest$iterations)),row.names=FALSE)
  invisible(summaries)
}
run_imputation_extensions <- function(args=commandArgs(trailingOnly=TRUE),script_dir=NULL) {
  if(is.null(script_dir)) {
    script <- grep("^--file=",commandArgs(FALSE),value=TRUE)
    if(length(script)!=1L) stop("Supply script_dir when called interactively.")
    script_dir <- dirname(normalizePath(sub("^--file=","",script)))
  }
  source(file.path(script_dir,"imputation_extensions.R"),local=environment(run_imputation_extensions))
  options <- mie_cli_options(args)
  selected <- if(options$scenario=="all")mie_scenarios()$scenario else options$scenario
  if(any(!selected %in% mie_scenarios()$scenario)) stop("Unknown scenario.")
  if(options$stage!="generate" && !dir.exists(options$output)) stop("Existing output root required.")
  dir.create(options$output,recursive=TRUE,showWarnings=FALSE)
  root <- normalizePath(options$output,winslash="/",mustWork=TRUE)
  records <- list(); results <- list()
  for(scenario in selected) {
    current <- options; current$scenario <- scenario; current$output <- file.path(root,"scenarios",scenario)
    # A batch invocation can inspect/fit scenarios at their independently extended counts.
    manifest_path <- file.path(current$output,"manifest.rds")
    if(options$stage!="generate" && file.exists(manifest_path)) {
      manifest <- readRDS(manifest_path); current$m <- manifest$m; current$iterations <- manifest$iterations
    }
    current <- current[!vapply(current,is.null,logical(1L))]
    current_args <- as.vector(rbind(paste0("--",names(current)),vapply(current,as.character,character(1L))))
    result <- tryCatch(list(value=run_imputation_scenario(current_args,script_dir),error=NULL),error=function(e)list(value=NULL,error=conditionMessage(e)))
    results[[scenario]] <- result$value
    records[[length(records)+1L]] <- data.frame(stage=options$stage,scenario=scenario,m=current$m,iterations=current$iterations,
      status=if(is.null(result$error))"completed" else "failed",reason=if(is.null(result$error))"" else result$error)
  }
  registry <- mie_bind(records)
  stamp <- format(Sys.time(),"%Y%m%dT%H%M%OS6"); stamp <- gsub("[^0-9T]","",stamp)
  write.csv(registry,file.path(root,paste0("stage_",options$stage,"_",stamp,".csv")),row.names=FALSE)
  if(any(registry$status=="failed")) stop("Some scenarios could not complete this stage; see the stage registry. Other scenarios were retained.")
  invisible(results)
}
if(sys.nframe()==0L) run_imputation_extensions()
