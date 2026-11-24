# Source this file from any working directory, then call run_experiment().
# Rscript path/to/airport_lasso.R runs the experiment directly.
.airport_model_dir <- local({
  # source() records the script path. With chdir=TRUE it already entered its folder.
  sources <- Filter(function(f) is.character(f$ofile) && length(f$ofile)==1L, sys.frames())
  if (length(sources)) {
    caller <- sources[[length(sources)]]
    script <- if (is.null(caller$owd)) caller$ofile else basename(caller$ofile)
  } else {
    args <- commandArgs(trailingOnly=FALSE)
    file_arg <- grep("^--file=", args, value=TRUE)
    script <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else args[match("-f", args)+1L]
    if (is.na(script)) script <- "airport_lasso.R"
    if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed=TRUE)
  }
  dirname(normalizePath(script, mustWork=TRUE))
})
if (!file.exists(file.path(.airport_model_dir,"preprocess_weather.R"))) {
  stop("Missing preprocess_weather.R in ", .airport_model_dir,
       ". Extract airport_lasso_code.zip and open common_design/airport_lasso.R.", call.=FALSE)
}
source(file.path(.airport_model_dir,"preprocess_weather.R"),local=TRUE)

source(file.path(.airport_model_dir,"score_mse.R"),local=TRUE)

load_solver <- function(directory) {
  src <- file.path(directory,"resources")
  library <- file.path(src,paste0("lasso",.Platform$dynlib.ext))
  if(!file.exists(library)) {
    old <- setwd(src)
    on.exit(setwd(old))
    if(system2(file.path(R.home("bin"),"R"),c("CMD","SHLIB","lasso.c"))!=0L)
      stop("Building the Lasso solver requires the R C toolchain.",call.=FALSE)
  }
  dyn.load(library)
}

forecast_airport <- function(model,five_rows,airport) {
  input <- prepare_forecast(model,five_rows,airport)
  forecasts <- lapply(model$fits[[airport]],function(f) {
    stopifnot(identical(colnames(input$raw),names(f$beta)))
    x <- apply_scaling(input$raw,f$scaling)
    pred <- input$normal[f$lead,f$target]+f$intercept+sum(x*f$beta)
    data.frame(airport,origin=input$origin,forecast_date=input$origin+f$lead,
               lead_days=f$lead,target=f$target,prediction_f=pred)
  })
  do.call(rbind,unname(forecasts))
}

fit_period <- function(end, eval_years, final=FALSE, project=load_weather_project(.airport_model_dir),
                       origins=forecast_origins(project$obs$dates,eval_years),remote_types="all") {
  dates <- project$obs$dates
  temp <- project$obs$temp
  stations <- project$obs$stations
  targets <- project$targets
  year_id <- as.integer(format(dates,"%Y"))
  workers <- if(.Platform$OS.type=="windows") 1L else as.integer(Sys.getenv("AIRPORT_WORKERS","4"))
  train <- which(year_id<=end)
  eligible <- forecast_origins(dates,unique(year_id[train]))
  stopifnot(length(eligible)>0L,length(origins)>0L,min(origins)>=max(train))
  daily <- prepare_days(project$obs,train)
  lambda <- 0.0625
  path_lambda <- lambda
  run_airport <- function(a) {
    raw <- airport_inputs(daily,a,remote_types)
    p <- ncol(raw)
    fits <- list()
    predictions <- list()
    training <- list()
    for (h in 1:5) {
      obs <- daily$anomaly[eligible+h,temperature_cols(a),drop=FALSE]
      valid <- lapply(1:3,function(j) which(is.finite(obs[,j])))
      for (idx in unique(valid)) {
        kinds <- which(vapply(valid,identical,logical(1),idx))
        rows <- eligible[idx]
        n <- length(rows)
        scaling <- learn_scaling(raw[rows,,drop=FALSE])
        x <- apply_scaling(raw[c(rows,origins),,drop=FALSE],scaling)
        xt <- x[seq_len(n),,drop=FALSE]
        newx <- x[n+seq_along(origins),,drop=FALSE]
        y <- obs[idx,kinds,drop=FALSE]
        means <- colMeans(y)
        yc <- sweep(y,2,means)
        g <- crossprod(xt)/n
        xy <- crossprod(xt,yc)/n
        for(i in seq_along(kinds)) {
          j <- kinds[i]
          solver <- .C("fit_path",g,as.double(xy[,i]),as.integer(p),as.double(path_lambda),
                       as.integer(length(path_lambda)),path=double(p*length(path_lambda)),
                       err=double(length(path_lambda)))
          stopifnot(max(solver$err)<1.01e-5)
          beta <- setNames(tail(solver$path,p),colnames(raw))
          pred <- daily$normal[origins+h,temperature_cols(a)[j]]+means[i]+as.vector(newx %*% beta)
          key <- paste(target_fields[j],h,sep="_")
          predictions[[key]] <- data.frame(airport=a,origin=dates[origins],lead_days=h,target=target_fields[j],
                     forecast_date=dates[origins+h],fit_end_year=end,
                     observed_f=temp[origins+h,temperature_cols(a)[j]],prediction_f=pred)
          if(final) {
            fits[[key]] <- list(lead=h,target=target_fields[j],intercept=unname(means[i]),
                                 beta=beta,scaling=scaling,training_n=n)
            training[[key]] <- data.frame(airport=a,lead_days=h,target=target_fields[j],
                             training_loss(yc[,i],as.vector(xt %*% beta)))
          }
        }
      }
    }
    list(pred=do.call(rbind,unname(predictions)),fits=fits,train=do.call(rbind,unname(training)))
  }
  all <- parallel::mclapply(targets,run_airport,mc.cores=workers,mc.preschedule=FALSE)
  stopifnot(all(vapply(all,function(v) is.list(v) && !is.null(v$pred),logical(1))))
  pred <- do.call(rbind,lapply(all,"[[","pred"))
  stopifnot(all(is.finite(pred$prediction_f)))
  model <- if(final) list(design="ordinary_lasso",remote_types=remote_types,lambda=lambda,
           training_start=min(dates[train]),training_end=max(dates[train]),stations=stations,
           season_start=dates[1],temp_season=daily$temp_season,pressure_scaling=daily$pressure_scaling,
           fits=setNames(lapply(all,"[[","fits"),targets))
  list(pred=pred,model=model,training=do.call(rbind,lapply(all,"[[","train")))
}

save_selected_predictors <- function(model,artifacts) {
  coefficients <- do.call(rbind,lapply(names(model$fits),function(a) {
    do.call(rbind,lapply(model$fits[[a]],function(f) {
      nz <- which(f$beta!=0)
      names_used <- names(f$beta)[nz]
      source <- sub("_.*$","",names_used)
      kind <- sub("_lag[0-9]+$","",sub("^[^_]+_","",names_used))
      lag <- as.integer(sub("^.*_lag","",names_used))
      data.frame(target_airport=rep(a,length(nz)),lead_days=rep(f$lead,length(nz)),
                 target=rep(f$target,length(nz)),source_airport=source,location=ifelse(source==a,"local","external"),
                 variable=kind,lag_days=lag,feature=names_used,
                 coefficient_per_training_sd=unname(f$beta[nz]))
    }))
  }))
  write.csv(coefficients,file.path(artifacts,"selected_predictors.csv"),row.names=FALSE)
}

run_experiment <- function(remote_types=c("all","mean"), directory=.airport_model_dir) {
  remote_types <- match.arg(remote_types)
  project <- load_weather_project(directory)
  load_solver(directory)
  artifacts <- file.path(result_directory(directory,remote_types),"artifacts")
  dir.create(artifacts,recursive=TRUE,showWarnings=FALSE)
  final <- fit_period(2012,2013:2014,TRUE,project,remote_types=remote_types)
  obj_path <- file.path(artifacts,"airport_lasso.rds")
  saveRDS(final$model,obj_path,version=3)
  model <- readRDS(obj_path)
  stopifnot(sum(lengths(model$fits))==225L)
  dates <- project$obs$dates
  # Check the five-row forecasting interface at the start, year boundary and end.
  evaluation_origins <- forecast_origins(dates,2013:2014)
  check <- unique(c(head(evaluation_origins,1),which(dates==as.Date("2013-12-31")),tail(evaluation_origins,1)))
  reconstructed <- do.call(rbind,lapply(check,function(z)
    do.call(rbind,lapply(project$targets,function(a) forecast_airport(model,project$raw[(z-4):z,,drop=FALSE],a)))))
  keys <- c("airport","origin","lead_days","target")
  paired <- merge(reconstructed,final$pred,by=keys)
  stopifnot(nrow(final$pred)==length(evaluation_origins)*225L,nrow(paired)==length(check)*225L,
            !anyDuplicated(paired[keys]),max(abs(paired$prediction_f.x-paired$prediction_f.y))<1e-8)
  save_selected_predictors(model,artifacts)
  # Assign all five leads to the fit available at the origin, including year-crossing windows.
  cv_origins <- forecast_origins(dates,2008:2012)
  issue_fold <- as.integer(format(dates[cv_origins+1L],"%Y"))
  folds <- lapply(2008:2012,function(year) {
    message("Scoring daily origins for CV fold ",year)
    fit_period(year-1,year,project=project,origins=cv_origins[issue_fold==year],remote_types=remote_types)$pred
  })
  run <- list(schedule="daily_5_to_5",remote_types=remote_types,training=final$training,cv=do.call(rbind,folds),evaluation=final$pred)
  pred_file <- file.path(artifacts,"forecast_errors_input.rds")
  saveRDS(run,pred_file)
  scores <- score_saved_predictions(pred_file)
  print(scores$summary)
  invisible(scores)
}

if (sys.nframe()==0L) run_experiment()
