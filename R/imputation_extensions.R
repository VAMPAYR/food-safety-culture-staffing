# Specification 1.4, priority 3. Definitions only: sourcing performs no analysis.
# No functions or globals from the historical imputation.R are replaced.
mie_variables <- function() c(paste0("Q", c(1,2,3,4,6,7,8,9,11,12,13,15,17,21,23,24,25,26,27)), "D1", "D2", "D5")
mie_seeds <- function() c(20261111L, 20261112L)
mie_scenarios <- function() data.frame(scenario = c("baseline", "both_lower", "both_higher", "busy_lower", "absence_lower"),
  delta_B = c(0,-1,1,-1,1)*log(2), delta_A = c(0,-1,1,1,-1)*log(2), stringsAsFactors = FALSE)
mie_bind <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
mie_with_rng <- function(state = NULL, code) {
  previous_kind <- RNGkind(); had_seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  if (had_seed) previous_seed <- get(".Random.seed", envir=.GlobalEnv)
  on.exit({do.call(RNGkind, as.list(previous_kind)); if (had_seed) assign(".Random.seed", previous_seed, envir=.GlobalEnv)
    else if (exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)) rm(".Random.seed", envir=.GlobalEnv)}, add=TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  if (!is.null(state)) assign(".Random.seed", state, envir=.GlobalEnv)
  force(code)
}
mie_streams <- function(seed, m) mie_with_rng(code={
  stopifnot(m >= 1L, m == as.integer(m)); set.seed(seed)
  streams <- vector("list",m); streams[[1L]] <- .Random.seed
  if (m>1L) for (i in 2:m) streams[[i]] <- parallel::nextRNGStream(streams[[i-1L]])
  streams
})
mie_encode <- function(data) {
  x <- data
  for (v in c("D1","D5")) if (is.factor(x[[v]])) x[[v]] <- as.integer(x[[v]])
  if (is.factor(x$D2)) x$D2 <- match(as.character(x$D2), c("false","true"))-1L
  x
}
mie_validate <- function(data) {
  required <- c("Observation","Site","Language",mie_variables())
  if (!is.data.frame(data) || !all(required %in% names(data))) stop("Missing enriched-PMM variables.")
  if (anyNA(data$Observation) || anyDuplicated(data$Observation)) stop("Invalid record identifiers.")
  site <- as.character(data$Site)
  if (anyNA(site) || !setequal(site,as.character(1:8))) stop("All eight observed sites are required.")
  if (anyNA(data$Language) || any(!as.character(data$Language) %in% c("English","Spanish"))) stop("Language must be observed and documented.")
  for (v in mie_variables()) {
    allowed <- if (v=="D2") 0:1 else 1:5
    if (!is.numeric(data[[v]]) || any(!is.na(data[[v]]) & !data[[v]] %in% allowed)) stop("Invalid decoded values: ",v)
    observed <- data[[v]][!is.na(data[[v]])]
    if (anyNA(data[[v]]) && length(observed)<5L) stop("Insufficient observed donor support: ",v)
  }
  invisible(TRUE)
}
mie_design <- function(completed, target) {
  if (!target %in% mie_variables()) stop("Unknown conditional target.")
  predictors <- setdiff(mie_variables(),target)
  columns <- list(`(Intercept)`=rep(1,nrow(completed)))
  for (v in predictors) for (level in if (v=="D2") 1L else 2:5)
    columns[[paste0(v,"_",level)]] <- as.numeric(completed[[v]]==level)
  for (site in 2:8) columns[[paste0("Site",site)]] <- as.numeric(as.character(completed$Site)==as.character(site))
  columns$LanguageSpanish <- as.numeric(as.character(completed$Language)=="Spanish")
  X <- do.call(cbind,columns); storage.mode(X) <- "double"
  if (any(!is.finite(X))) stop("Incomplete or invalid conditional predictors.")
  X
}
mie_posterior <- function(X,z) {
  if (!is.matrix(X) || nrow(X)!=length(z) || !length(z) || any(!is.finite(c(X,z)))) stop("Invalid working-regression input.")
  if (any(z<0 | z>1)) stop("Conditional target must be on [0,1].")
  q <- ncol(X); prior_mean <- c(.5,rep(0,q-1L)); gram <- crossprod(X)
  precision <- gram + diag(.01,q); R <- chol(precision)
  mean <- drop(backsolve(R,forwardsolve(t(R),crossprod(X,z)+.01*prior_mean)))
  shape <- 2+length(z)/2
  scale <- .04+.5*(sum((z-drop(X%*%mean))^2)+.01*sum((mean-prior_mean)^2))
  sigma2 <- scale/rgamma(1L,shape=shape,rate=1)
  draw <- mean+sqrt(sigma2)*drop(backsolve(R,rnorm(q)))
  if (any(!is.finite(c(mean,scale,sigma2,draw))) || scale<=0 || sigma2<=0) stop("Invalid conditional posterior draw.")
  list(mean=mean,draw=draw,sigma2=sigma2,precision=precision,shape=shape,scale=scale,chol=R,gram=gram)
}
mie_weights <- function(values,delta=0) {
  if (!length(values) || any(!is.finite(c(values,delta))) || length(delta)!=1L) stop("Invalid donor tilt.")
  log_weight <- delta*(values-3); w <- exp(log_weight-max(log_weight)); w/sum(w)
}
mie_draw_donors <- function(pred_obs,pred_mis,original_values,delta=0) {
  if (length(pred_obs)<5L || length(original_values)!=length(pred_obs) || any(!is.finite(c(pred_obs,pred_mis,original_values)))) stop("Invalid five-donor matching input.")
  selected <- integer(length(pred_mis)); details <- vector("list",length(pred_mis))
  sets <- matrix(integer(length(pred_mis)*5L),ncol=5L)
  for (i in seq_along(pred_mis)) {
    index <- order(abs(pred_obs-pred_mis[i]),runif(length(pred_obs)))[1:5]
    values <- original_values[index]; w <- mie_weights(values,delta)
    # A single uniform draw is also used when delta=0, enabling common streams.
    position <- which(cumsum(w)>=runif(1L))[1L]
    if (is.na(position)) position <- 5L
    selected[i] <- index[position]; sets[i,] <- index
    counts <- tabulate(match(values,0:5),nbins=6L)
    details[[i]] <- data.frame(recipient=i,unique_categories=length(unique(values)),single_category=length(unique(values))==1L,
      conditional_mean_shift=sum(values*w)-mean(values),agreement_pp_shift=100*(sum(w[values>=4])-mean(values>=4)),
      candidate_mean=mean(values),tilted_mean=sum(values*w),n0=counts[1],n1=counts[2],n2=counts[3],n3=counts[4],n4=counts[5],n5=counts[6])
  }
  list(selected=selected,details=mie_bind(details),sets=sets)
}
mie_update <- function(completed,original,target,scenario,terminal=FALSE) {
  observed <- which(!is.na(original[[target]])); missing <- which(is.na(original[[target]]))
  X <- mie_design(completed,target)
  z <- original[[target]][observed]; if (target!="D2") z <- (z-1)/4
  post <- mie_posterior(X[observed,,drop=FALSE],z)
  delta <- if (target=="Q11") scenario$delta_B else if (target=="Q12") scenario$delta_A else 0
  donors <- mie_draw_donors(drop(X[observed,,drop=FALSE]%*%post$mean),drop(X[missing,,drop=FALSE]%*%post$draw),original[[target]][observed],delta)
  selected <- observed[donors$selected]; completed[[target]][missing] <- original[[target]][selected]
  donors$details$record <- as.character(original$Observation[missing])
  donors$details$selected_record <- as.character(original$Observation[selected])
  donors$details$target <- target
  if (terminal) for (k in 1:5) donors$details[[paste0("donor",k)]] <- as.character(original$Observation[observed[donors$sets[,k]]])
  diagnostic <- NULL
  if (terminal) {
    eig <- eigen(post$gram,symmetric=TRUE,only.values=TRUE)$values
    diagnostic <- data.frame(target=target,n_observed=length(observed),columns=ncol(X),rank=qr(X[observed,,drop=FALSE])$rank,
      minimum_data_eigenvalue=min(eig),maximum_data_eigenvalue=max(eig),precision_condition=(max(eig)+.01)/(min(eig)+.01),
      maximum_prior_precision_fraction=.01/(min(eig)+.01),posterior_shape=post$shape,posterior_scale=post$scale)
  }
  list(data=completed,selected=selected,details=donors$details,diagnostic=diagnostic)
}
mie_initialize <- function(original,stream,scenario,seed,chain) mie_with_rng(stream,{
  completed <- original
  incomplete <- mie_variables()[vapply(original[mie_variables()],anyNA,logical(1L))]
  for (v in incomplete) {missing <- is.na(original[[v]]); completed[[v]][missing] <- sample(original[[v]][!missing],sum(missing),replace=TRUE)}
  list(version="1.4",scenario=scenario,seed=as.integer(seed),chain=as.integer(chain),iteration=0L,status="initialized",
    data=completed,stream=.Random.seed,trace=list(),donor_counts=setNames(lapply(incomplete,function(v) integer(nrow(original))),incomplete),
    staffing_updates=list(),terminal_donors=list(),terminal_regression=list(),failure=NULL)
})
mie_save <- function(object,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temporary <- paste0(path,".tmp"); saveRDS(object,temporary,version=3L,compress=FALSE)
  # Keep a recoverable previous checkpoint while installing the new file.
  previous <- paste0(path,".previous")
  if (file.exists(previous)) unlink(previous)
  if (file.exists(path) && !file.rename(path,previous)) stop("Cannot preserve previous checkpoint: ",path)
  if (!file.rename(temporary,path)) stop("Cannot install checkpoint; previous remains: ",path)
  if (file.exists(previous)) unlink(previous)
  invisible(path)
}
mie_advance <- function(state,original,iterations,path=NULL) {
  if (state$status=="failed") return(state)
  if (state$iteration>iterations) stop("Cannot shorten a chain.")
  if (state$iteration==iterations) return(state)
  mie_with_rng(state$stream,{
    incomplete <- names(state$donor_counts)
    for (iteration in seq.int(state$iteration+1L,iterations)) {
      terminal <- iteration==iterations; working <- state$data
      target <- NA_character_
      attempt <- tryCatch({
        for (target in incomplete) {
          update <- mie_update(working,original,target,state$scenario,terminal)
          working <- update$data
          state$donor_counts[[target]] <- state$donor_counts[[target]]+tabulate(update$selected,nrow(original))
          values <- working[[target]][is.na(original[[target]])]; counts <- tabulate(match(values,0:5),nbins=6)
          state$trace[[length(state$trace)+1L]] <- data.frame(iteration=iteration,target=target,n_imputed=length(values),mean=mean(values),
            variance=if(length(values)>1) var(values) else NA_real_,n0=counts[1],n1=counts[2],n2=counts[3],n3=counts[4],n4=counts[5],n5=counts[6])
          if (target %in% c("Q11","Q12")) {row <- update$details[,!grepl("^donor[1-5]$",names(update$details)),drop=FALSE]; row$iteration <- iteration; state$staffing_updates[[length(state$staffing_updates)+1L]] <- row}
          if (terminal) {state$terminal_donors[[target]] <- update$details; state$terminal_regression[[target]] <- update$diagnostic}
        }
        NULL
      },error=function(e) conditionMessage(e))
      state$stream <- .Random.seed
      if (!is.null(attempt)) {
        state$status <- "failed"; state$failure <- list(iteration=iteration,target=target,reason=attempt,partial_data=working)
        if (!is.null(path)) mie_save(state,path)
        return(state)
      }
      state$data <- working; state$iteration <- iteration; state$status <- if(terminal) "complete" else "running"
      if (!is.null(path)) mie_save(state,path)
    }
    state
  })
}
mie_analysis_data <- function(original,completed,outcome) {
  if (!outcome %in% c(4L,6L)) stop("Unsupported outcome.")
  completed[!is.na(original[[paste0("Q",outcome)]]),,drop=FALSE]
}
mie_fit_one <- function(original,state,outcome) {
  if (!identical(state$status,"complete")) stop("Completed chain required.")
  completed <- mie_analysis_data(original,state$data,outcome)
  prepared <- derive_scores(completed,outcome=outcome,management="stored")
  id <- sprintf("enriched_%s_seed%d_chain%d_iter%d_Q%d",state$scenario$scenario,state$seed,state$chain,state$iteration,outcome)
  fit <- tryCatch(fit_ordinal(prepared,"Y",c("L_c","E_c","R_c","M_c","B_c","A_c"),site="Site",model_id=id),
    error=function(e) list(valid=FALSE,status="error",reason=conditionMessage(e)))
  coefficients <- if (!is.null(fit$coef)) fit$coef else fit$coefficients
  valid <- isTRUE(fit$valid) && all(c("B_c","A_c") %in% names(coefficients)) && all(is.finite(coefficients[c("B_c","A_c")]))
  values <- if (valid) c(beta_B=unname(coefficients["B_c"]),beta_A=unname(coefficients["A_c"])) else c(beta_B=NA_real_,beta_A=NA_real_)
  if (outcome==6L) values <- c(values,d=unname(values["beta_A"]-values["beta_B"]))
  list(scenario=state$scenario$scenario,seed=state$seed,chain=state$chain,iteration=state$iteration,outcome=outcome,n=nrow(completed),
    valid=valid,reason=if(valid) "" else paste(c(fit$reason,fit$status),collapse="; "),values=values,fit=fit)
}
mie_summarize_fits <- function(registry,m,scenarios=mie_scenarios()$scenario,iterations=NULL) {
  estimates <- covariance <- availability <- checks <- list()
  for (scenario in scenarios) {
    scenario_iterations <- unique(vapply(Filter(function(x)identical(x$scenario,scenario),registry),function(x)as.integer(x$iteration),integer(1L)))
    intended_iteration <- if(!is.null(iterations))as.integer(iterations) else if(length(scenario_iterations)==1L)scenario_iterations else NA_integer_
    scenario_estimates <- scenario_covariance <- list(); outcome_available <- setNames(rep(TRUE,2),c("6","4"))
    for (seed in mie_seeds()) for (outcome in c(6L,4L)) {
      rows <- Filter(function(x) identical(x$scenario,scenario) && x$seed==seed && x$outcome==outcome,registry)
      record_valid <- vapply(rows,function(x) {
        expected <- if(outcome==6L)c("beta_B","beta_A","d") else c("beta_B","beta_A")
        isTRUE(x$valid) && !is.na(intended_iteration) && x$iteration==intended_iteration &&
          identical(names(x$values),expected) && all(is.finite(x$values)) &&
          (outcome!=6L || abs(x$values["d"]-(x$values["beta_A"]-x$values["beta_B"]))<1e-12)
      },logical(1L))
      valid <- length(rows)==m && setequal(vapply(rows,function(x)x$chain,integer(1L)),seq_len(m)) &&
        !anyDuplicated(vapply(rows,function(x)x$chain,integer(1L))) && all(record_valid)
      n_valid <- sum(record_valid)
      availability[[length(availability)+1L]] <- data.frame(scenario=scenario,seed=seed,outcome=outcome,required=m,attempted=length(rows),valid=n_valid,failed=length(rows)-n_valid,status=if(valid)"available" else "unavailable")
      outcome_available[as.character(outcome)] <- outcome_available[as.character(outcome)] && valid
      if (!valid) next
      values <- do.call(rbind,lapply(rows,function(x)x$values)); average <- colMeans(values); between <- cov(values)
      scenario_estimates[[length(scenario_estimates)+1L]] <- data.frame(scenario=scenario,seed=seed,outcome=outcome,component=names(average),m=m,
        estimate=unname(average),exponentiated=exp(unname(average)),exponentiated_quantity=ifelse(names(average)=="d","ratio_of_odds_ratios","odds_ratio"),
        between_imputation_variance=unname(diag(between)),monte_carlo_error=sqrt(unname(diag(between))/m))
      index <- expand.grid(component_row=colnames(between),component_column=colnames(between),stringsAsFactors=FALSE)
      scenario_covariance[[length(scenario_covariance)+1L]] <- data.frame(scenario=scenario,seed=seed,outcome=outcome,index,between_imputation_covariance=as.vector(between))
    }
    e <- mie_bind(scenario_estimates); cv <- mie_bind(scenario_covariance)
    if (nrow(e)) {e <- e[outcome_available[as.character(e$outcome)],,drop=FALSE]; estimates[[length(estimates)+1L]] <- e}
    if (nrow(cv)) covariance[[length(covariance)+1L]] <- cv[outcome_available[as.character(cv$outcome)],,drop=FALSE]
    if (nrow(e)) for (outcome in c(6L,4L)[outcome_available]) {
      first <- e[e$seed==mie_seeds()[1]&e$outcome==outcome,,drop=FALSE]
      second <- e[e$seed==mie_seeds()[2]&e$outcome==outcome,,drop=FALSE]; second <- second[match(first$component,second$component),,drop=FALSE]
      difference <- abs(first$estimate-second$estimate); mc <- pmax(first$monte_carlo_error,second$monte_carlo_error)
      checks[[length(checks)+1L]] <- data.frame(scenario=scenario,outcome=outcome,component=first$component,absolute_seed_difference=difference,
        maximum_monte_carlo_error=mc,exceeds_tolerance=difference>.02|mc>.01)
    }
  }
  checks <- mie_bind(checks)
  list(estimates=mie_bind(estimates),covariance=mie_bind(covariance),availability=mie_bind(availability),monte_carlo_checks=checks,
    extension_required=nrow(checks)>0 && any(checks$exceeds_tolerance),at_chain_cap=m==200L)
}

mie_chain_path <- function(outdir,scenario,seed,chain) file.path(outdir,"chains",scenario,paste0("seed",seed),sprintf("chain%04d.rds",chain))
mie_tasks <- function(outdir,m,iterations,scenarios=mie_scenarios()$scenario) {
  tasks <- expand.grid(scenario=scenarios,seed=mie_seeds(),chain=seq_len(m),stringsAsFactors=FALSE)
  tasks$iterations <- iterations
  tasks$path <- mapply(mie_chain_path,MoreArgs=list(outdir=outdir),scenario=tasks$scenario,seed=tasks$seed,chain=tasks$chain,USE.NAMES=FALSE)
  tasks
}
mie_generation_worker <- function(task) {
  path <- task$path
  tryCatch({
    state <- if (file.exists(path)) readRDS(path) else mie_initialize(mie_worker_original,mie_worker_streams[[as.character(task$seed)]][[task$chain]],
      mie_scenarios()[match(task$scenario,mie_scenarios()$scenario),,drop=FALSE],task$seed,task$chain)
    if (state$scenario$scenario!=task$scenario || state$seed!=task$seed || state$chain!=task$chain) stop("Chain identity mismatch.")
    state <- mie_advance(state,mie_worker_original,task$iterations,path)
    data.frame(task,status=state$status,iteration=state$iteration,reason=if(is.null(state$failure))"" else state$failure$reason,check.names=FALSE)
  },error=function(e) {
    failure <- list(task=task,reason=conditionMessage(e),time=as.character(Sys.time()))
    mie_save(failure,paste0(path,".failure.rds"))
    data.frame(task,status="failed",iteration=NA_integer_,reason=conditionMessage(e),check.names=FALSE)
  })
}
mie_parallel <- function(tasks,fun,script_dir,original,streams=NULL,workers=4L) {
  if (!workers %in% 1:4) stop("Use between one and four PSOCK workers.")
  cluster <- parallel::makePSOCKcluster(workers,outfile="")
  on.exit(parallel::stopCluster(cluster),add=TRUE)
  parallel::clusterCall(cluster,function(directory,original,streams) {
    Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    source(file.path(directory,"data_prep.R"),local=.GlobalEnv)
    source(file.path(directory,"model_engine.R"),local=.GlobalEnv)
    source(file.path(directory,"imputation_extensions.R"),local=.GlobalEnv)
    assign("mie_worker_original",original,envir=.GlobalEnv); assign("mie_worker_streams",streams,envir=.GlobalEnv)
    NULL
  },script_dir,original,streams)
  parallel::parLapplyLB(cluster,lapply(seq_len(nrow(tasks)),function(i)as.list(tasks[i,,drop=FALSE])),fun)
}
mie_write_diagnostics <- function(outdir,original,m,iterations,scenarios=mie_scenarios()$scenario) {
  tasks <- mie_tasks(outdir,m,iterations,scenarios)
  tables <- setNames(vector("list",10),c("chains","trace","distributions","donor_reuse","terminal_donors","terminal_regression","staffing_updates","associations","staffing_completed","analysis_missing_staffing"))
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i,,drop=FALSE]
    state <- if(file.exists(task$path)) tryCatch(readRDS(task$path),error=function(e)NULL) else NULL
    prefix <- task[,c("scenario","seed","chain"),drop=FALSE]
    identity_ok <- !is.null(state) && identical(state$scenario$scenario,task$scenario) && state$seed==task$seed && state$chain==task$chain
    sidecar_failed <- file.exists(paste0(task$path,".failure.rds"))
    complete <- identity_ok && identical(state$status,"complete") && state$iteration==iterations && !sidecar_failed
    tables$chains[[i]] <- data.frame(prefix,required_iterations=iterations,iteration=if(is.null(state))NA_integer_ else state$iteration,
      status=if(is.null(state))"missing_or_unreadable" else state$status,complete=complete,
      reason=if(!identity_ok)"missing_or_mismatched_checkpoint" else if(sidecar_failed)"retained_generation_failure" else if(is.null(state$failure))"" else state$failure$reason)
    if(!identity_ok) next
    append_table <- function(name,x) {if(nrow(x)) tables[[name]][[length(tables[[name]])+1L]] <<- cbind(prefix[rep(1L,nrow(x)),,drop=FALSE],x)}
    append_table("trace",mie_bind(state$trace)); append_table("terminal_donors",mie_bind(state$terminal_donors))
    append_table("terminal_regression",mie_bind(state$terminal_regression)); append_table("staffing_updates",mie_bind(state$staffing_updates))
    for(v in names(state$donor_counts)) append_table("donor_reuse",data.frame(target=v,record=as.character(original$Observation),times_selected=state$donor_counts[[v]]))
    if(!complete) next
    for(v in mie_variables()) for(site in c("all",as.character(1:8))) for(kind in c("observed","imputed")) {
      rows <- (site=="all" | as.character(original$Site)==site) & (is.na(original[[v]])==(kind=="imputed"))
      values <- state$data[[v]][rows]; support <- if(v=="D2")0:1 else 1:5
      append_table("distributions",data.frame(target=v,site=site,kind=kind,category=support,n=tabulate(match(values,support),length(support)),denominator=length(values)))
    }
    for(a in seq_along(mie_variables())) if(a<length(mie_variables())) for(b in seq.int(a+1L,length(mie_variables()))) {
      va <- mie_variables()[a]; vb <- mie_variables()[b]
      for(scope in c("all","at_least_one_imputed")) {
        rows <- if(scope=="all")rep(TRUE,nrow(original)) else is.na(original[[va]])|is.na(original[[vb]])
        usable <- sum(rows)>1 && sd(state$data[[va]][rows])>0 && sd(state$data[[vb]][rows])>0
        append_table("associations",data.frame(item_a=va,item_b=vb,scope=scope,n=sum(rows),correlation=if(usable)cor(state$data[[va]][rows],state$data[[vb]][rows]) else NA_real_))
      }
    }
    for(outcome in c(6L,4L)) for(v in c("Q11","Q12")) {
      rows <- !is.na(original[[paste0("Q",outcome)]])
      missing <- rows & is.na(original[[v]])
      values <- state$data[[v]][rows]
      append_table("staffing_completed",data.frame(outcome=outcome,item=v,n=sum(rows),category=1:5,count=tabulate(values,5),
        missing_cells=sum(missing),imputed_mean=if(any(missing))mean(state$data[[v]][missing])else NA_real_))
    }
  }
  for(outcome in c(6L,4L)) for(v in c("Q11","Q12")) tables$analysis_missing_staffing[[length(tables$analysis_missing_staffing)+1L]] <-
    data.frame(outcome=outcome,item=v,n_outcome_observed=sum(!is.na(original[[paste0("Q",outcome)]])),missing_staffing=sum(!is.na(original[[paste0("Q",outcome)]])&is.na(original[[v]])))
  directory <- file.path(outdir,sprintf("diagnostics_m%d_i%d",m,iterations)); dir.create(directory,showWarnings=FALSE)
  tables <- lapply(tables,mie_bind)
  for(name in names(tables)) write.csv(tables[[name]],file.path(directory,paste0(name,".csv")),row.names=FALSE)
  staffing <- tables$staffing_completed
  if(nrow(staffing)) {
    means <- unique(staffing[,c("scenario","seed","chain","outcome","item","imputed_mean")])
    baseline <- means[means$scenario=="baseline",,drop=FALSE]
    baseline_file <- NULL
    if(!nrow(baseline)) {
      baseline_directory <- file.path(dirname(outdir),"baseline")
      baseline_manifest <- file.path(baseline_directory,"manifest.rds")
      if(file.exists(baseline_manifest)) {
        bm <- readRDS(baseline_manifest)
        baseline_file <- file.path(baseline_directory,sprintf("diagnostics_m%d_i%d",bm$m,bm$iterations),"staffing_completed.csv")
        if(file.exists(baseline_file)) baseline <- unique(read.csv(baseline_file,stringsAsFactors=FALSE)[,c("scenario","seed","chain","outcome","item","imputed_mean")])
      }
    }
    key <- function(d)paste(d$seed,d$outcome,d$item,sep="|")
    baseline_average <- if(nrow(baseline))tapply(baseline$imputed_mean,key(baseline),mean) else numeric()
    means$baseline_imputed_mean <- unname(baseline_average[key(means)])
    means$baseline_reference_md5 <- if(!is.null(baseline_file)&&file.exists(baseline_file))unname(tools::md5sum(baseline_file)) else if(nrow(baseline))"same_scenario" else NA_character_
    means$realized_mean_difference <- means$imputed_mean-means$baseline_imputed_mean
    write.csv(means,file.path(directory,"staffing_realized_differences.csv"),row.names=FALSE)
    by_seed <- aggregate(cbind(imputed_mean,baseline_imputed_mean,realized_mean_difference)~scenario+seed+outcome+item,data=means,FUN=mean,na.action=na.pass)
    write.csv(by_seed,file.path(directory,"staffing_realized_mean_by_seed.csv"),row.names=FALSE)
  }
  files <- list.files(directory,pattern="\\.csv$",full.names=TRUE)
  files <- files[!basename(files) %in% c("diagnostic_md5.csv","chain_md5.csv")]
  diagnostic_hashes <- data.frame(file=basename(files),md5=unname(tools::md5sum(files)))
  write.csv(diagnostic_hashes,file.path(directory,"diagnostic_md5.csv"),row.names=FALSE)
  tasks$md5 <- ifelse(file.exists(tasks$path),unname(tools::md5sum(tasks$path)),NA_character_)
  write.csv(tasks,file.path(directory,"chain_md5.csv"),row.names=FALSE)
  fingerprint_files <- file.path(directory,c("diagnostic_md5.csv","chain_md5.csv"))
  list(directory=directory,complete=all(tables$chains$complete),fingerprint=unname(tools::md5sum(fingerprint_files)),files=fingerprint_files)
}
mie_load_saved_fit <- function(task) {
  failure_path <- paste0(task$fit_path,".failure.rds")
  path <- if(file.exists(failure_path))failure_path else task$fit_path
  if(!file.exists(path)) return(NULL)
  fit <- readRDS(path)
  if(!identical(fit$scenario,task$scenario) || fit$seed!=task$seed || fit$chain!=task$chain || fit$iteration!=task$iterations || fit$outcome!=task$outcome) stop("Cached fit identity mismatch.")
  fit
}
mie_validate_diagnostics <- function(diagnostics) {
  if(!identical(unname(tools::md5sum(diagnostics$files)),diagnostics$fingerprint)) stop("Diagnostic index has changed.")
  chains <- read.csv(file.path(diagnostics$directory,"chain_md5.csv"),stringsAsFactors=FALSE)
  if(anyNA(chains$md5)||!identical(unname(tools::md5sum(chains$path)),chains$md5)) stop("Chains changed since diagnostic review.")
  content <- read.csv(file.path(diagnostics$directory,"diagnostic_md5.csv"),stringsAsFactors=FALSE)
  if(!identical(unname(tools::md5sum(file.path(diagnostics$directory,content$file))),content$md5)) stop("Diagnostic content changed.")
  invisible(TRUE)
}
mie_require_review <- function(manifest_path,diagnostics,review_path,require_accept=TRUE) {
  if(!file.exists(review_path)) stop("A recorded diagnostic review is required before outcome fitting.")
  review <- readRDS(review_path); mie_validate_diagnostics(diagnostics)
  if(!identical(review$manifest_md5,unname(tools::md5sum(manifest_path))) || !identical(review$diagnostic_fingerprint,diagnostics$fingerprint)) stop("Review does not bind the current generation.")
  if(require_accept && (review$decision!="accept" || !isTRUE(diagnostics$complete))) stop("Current generation has no valid accepted review.")
  invisible(review)
}
mie_fit_worker <- function(task) {
  failure_path <- paste0(task$fit_path,".failure.rds")
  tryCatch({
    cached <- mie_load_saved_fit(task); if(!is.null(cached))return(cached)
    state <- readRDS(task$path)
    if(state$scenario$scenario!=task$scenario || state$seed!=task$seed || state$chain!=task$chain || state$iteration!=task$iterations) stop("Fit source checkpoint identity mismatch.")
    fit <- mie_fit_one(mie_worker_original,state,task$outcome)
    mie_save(fit,task$fit_path)
    fit
  },error=function(e) {
    values <- c(beta_B=NA_real_,beta_A=NA_real_); if(task$outcome==6L)values <- c(values,d=NA_real_)
    failure <- list(scenario=task$scenario,seed=task$seed,chain=task$chain,iteration=task$iterations,outcome=task$outcome,n=NA_integer_,
      valid=FALSE,reason=conditionMessage(e),values=values,fit=NULL)
    if(file.exists(failure_path)) {
      backup <- paste0(failure_path,".rejected_",gsub("[^0-9]","",format(Sys.time(),"%Y%m%d%H%M%OS6")))
      if(!file.copy(failure_path,backup,overwrite=FALSE)) stop("Could not preserve rejected failure artifact.")
    }
    mie_save(failure,failure_path); failure
  })
}
mie_fit_tasks <- function(outdir,m,iterations,scenarios=mie_scenarios()$scenario) {
  tasks <- mie_tasks(outdir,m,iterations,scenarios); tasks <- tasks[rep(seq_len(nrow(tasks)),each=2L),,drop=FALSE]
  tasks$outcome <- rep(c(6L,4L),nrow(tasks)/2L)
  tasks$fit_path <- file.path(outdir,"fits",sprintf("iteration%04d",iterations),
    sprintf("%s_seed%d_chain%04d_Q%d.rds",tasks$scenario,tasks$seed,tasks$chain,tasks$outcome))
  tasks
}
mie_write_fit_summary <- function(registry,outdir,m,iterations,scenarios=mie_scenarios()$scenario) {
  directory <- file.path(outdir,sprintf("summaries_m%d_i%d",m,iterations)); dir.create(directory,showWarnings=FALSE)
  summaries <- mie_summarize_fits(registry,m,scenarios,iterations=iterations)
  for(name in c("estimates","covariance","availability","monte_carlo_checks")) write.csv(summaries[[name]],file.path(directory,paste0(name,".csv")),row.names=FALSE)
  flat <- lapply(registry,function(x)data.frame(scenario=x$scenario,seed=x$seed,chain=x$chain,iteration=x$iteration,outcome=x$outcome,n=x$n,valid=x$valid,reason=x$reason))
  write.csv(mie_bind(flat),file.path(directory,"model_registry.csv"),row.names=FALSE)
  mie_save(summaries,file.path(directory,"summary.rds")); summaries
}
