.comparison_dir <- local({
  # source() records the script path. With chdir=TRUE it already entered its folder.
  sources <- Filter(function(f) is.character(f$ofile) && length(f$ofile)==1L, sys.frames())
  if (length(sources)) {
    caller <- sources[[length(sources)]]
    script <- if (is.null(caller$owd)) caller$ofile else basename(caller$ofile)
  } else {
    args <- commandArgs(trailingOnly=FALSE)
    file_arg <- grep("^--file=", args, value=TRUE)
    script <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else args[match("-f", args)+1L]
    if (is.na(script)) script <- "summarize_results.R"
    if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed=TRUE)
  }
  dirname(normalizePath(script, mustWork=TRUE))
})
source(file.path(.comparison_dir,"airport_lasso.R"),local=TRUE)

# Recompute the report tables from saved forecast errors and split-level losses.
summarize_results <- function(directory=.comparison_dir) {
  root <- file.path(dirname(directory),"main_results")
  dest <- file.path(root,"experiments")
  write <- function(x,name) write.csv(x,file.path(dest,paste0(name,".csv")),row.names=FALSE)
  loss <- read.csv(file.path(dest,"loss_by_fit_fold_and_penalty.csv"))
  cv <- subset(loss,split=="cv" & method %in% c("ridge","individual"))
  tuning <- aggregate(cbind(n,sse)~method+lambda,cv,sum)
  tuning$cv_mse <- tuning$sse/tuning$n
  tuning$selected <- with(tuning,ave(cv_mse,method,FUN=function(x)x==min(x)))==1
  write(tuning,"lambda_cv")
  labels <- subset(tuning,selected,c(method,lambda))
  chosen <- merge(labels,loss,by=c("method","lambda"))
  chosen$method <- paste0(chosen$method,"_tuned")
  rows <- rbind(subset(loss,!method %in% c("individual","ridge")),chosen)
  cols <- c("method","airport","lead_days","target","fit_end","split","n","sse")
  project <- load_weather_project(directory)
  dates <- project$obs$dates
  train <- which(as.integer(format(dates,"%Y"))<=2012)
  eligible <- forecast_origins(dates,2005:2012)
  monthly <- eligible[format(dates[eligible],"%d")=="15"]
  model <- readRDS(file.path(root,"baseline","artifacts","airport_lasso.rds"))
  daily <- prepare_days(project$obs,fitted=model)
  read_run <- function(variant) {
    run <- readRDS(file.path(root,variant,"artifacts","forecast_errors_input.rds"))
    fold_loss <- function(z) {
      z$n <- as.integer(is.finite(z$observed_f))
      z$sse <- ifelse(z$n==1,(z$observed_f-z$prediction_f)^2,0)
      aggregate(cbind(n,sse)~airport+lead_days+target+fit_end_year+split,z,sum) |>
        setNames(c("airport","lead_days","target","fit_end","split","n","sse"))
    }
    z <- rbind(transform(run$cv,split="cv"),transform(run$evaluation,split=format(forecast_date,"%Y")))
    training <- transform(run$training,split="training",fit_end=2012)
    one <- rbind(training,fold_loss(z))
    one$method <- variant
    if(variant=="baseline") {
      tr <- do.call(rbind,lapply(names(model$fits),function(a) {
        raw <- airport_inputs(daily,a)[monthly,,drop=FALSE]
        do.call(rbind,lapply(model$fits[[a]],function(f) {
          pred <- as.vector(apply_scaling(raw,f$scaling) %*% f$beta)+f$intercept
          gold <- daily$anomaly[monthly+f$lead,paste0(a,"_",f$target,".TemperatureF")]
          ok <- is.finite(gold)
          data.frame(airport=a,lead_days=f$lead,target=f$target,training_loss(gold[ok],pred[ok]))
        }))
      }))
      m <- rbind(transform(tr,split="training",fit_end=2012),fold_loss(z[format(z$origin,"%d")=="15",]))
      m$method <- "baseline_day15"
      one <- rbind(one,m)
    }
    one[,cols]
  }
  rows <- rbind(rows[,cols],read_run("baseline"),read_run("mean_only"))
  tbl <- function(data,keys) {
    test <- data[data$split %in% c("2013","2014"),]
    test$split <- "pooled"
    tmp <- rbind(data,test)
    totals <- aggregate(tmp[c("n","sse")],tmp[c(keys,"split")],sum)
    totals$mse <- totals$sse/totals$n
    wide <- reshape(totals[c(keys,"split","mse","n")],direction="wide",idvar=keys,timevar="split")
    wide[do.call(order,unname(as.list(wide[keys]))),c(keys,unlist(lapply(c("training","cv","2013","2014","pooled"),function(k) paste(c("mse","n"),k,sep="."))))]
  }
  write(tbl(rows,"method"),"model_scores")
  write(tbl(rows,c("method","lead_days","target")),"by_lead_target")
  write(tbl(rows,c("method","lead_days")),"by_lead")
  write(tbl(rows,c("method","target")),"by_target")
  splits <- subset(rows,split=="cv")
  fold <- aggregate(cbind(n,sse)~method+fit_end,splits,sum)
  fold$cv_mse <- fold$sse/fold$n
  write(fold,"validation_folds")
  key <- c("mean_only","baseline")
  pair <- function(keys) {
    t <- tbl(subset(rows,method %in% key),c("method",keys))
    groups <- if(length(keys)) unique(t[keys]) else data.frame(dummy=1)
    get <- function(k) t[t$method==k,,drop=FALSE]
    base <- get(key[1]);full <- get(key[2])
    stopifnot(identical(unname(as.matrix(base[keys])),unname(as.matrix(full[keys]))))
    for(sp in c("training","cv","2013","2014","pooled")) {
      value <- paste0("mse.",sp)
      stopifnot(all(base[[paste0("n.",sp)]]==full[[paste0("n.",sp)]]))
      groups[[paste0(sp,"_mean_only")]] <- base[[value]]
      groups[[paste0(sp,"_all_three")]] <- full[[value]]
      groups[[paste0(sp,"_difference")]] <- full[[value]]-base[[value]]
    }
    if(!length(keys)) groups$dummy <- NULL
    groups
  }
  write(pair(character()),"external_temperature_comparison")
  write(pair("target"),"external_temperature_by_target")
  write(pair("lead_days"),"external_temperature_by_lead")
  write(pair(c("lead_days","target")),"external_temperature_by_lead_target")
  selection <- function(variant) {
    sel <- read.csv(file.path(root,variant,"artifacts","selected_predictors.csv"))
    coord <- project$coords
    sel$distance_km <- NA_real_
    for(a in unique(sel$target_airport)) {
      matched <- which(sel$target_airport==a)
      sel$distance_km[matched] <- airport_distances(coord,a)[match(sel$source_airport[matched],coord$code)]
    }
    groups <- split(sel,interaction(sel$target_airport,sel$target,sel$lead_days))
    per_fit <- do.call(rbind,lapply(groups,function(z) {
      external <- z[z$location=="external",]
      distances <- external$distance_km[!duplicated(external$source_airport)]
      data.frame(airport=z$target_airport[1],target=z$target[1],lead_days=z$lead_days[1],
                 kept=nrow(z),local_kept=sum(z$location=="local"),external_airports=length(distances),
                 mapped_airports=sum(is.finite(distances)),median_distance_km=median(distances,na.rm=TRUE),
                 lag4_external_mean=sum(external$variable=="Tmean_anomaly" & external$lag_days==4),
                 lag4_external_temp=sum(grepl("^T",external$variable) & external$lag_days==4))
    }))
    stopifnot(nrow(per_fit)==225L)
    write(per_fit,paste0(variant,"_selection_by_fit"))
    measures <- c("kept","local_kept","external_airports","mapped_airports","median_distance_km","lag4_external_mean","lag4_external_temp")
    write(aggregate(per_fit[measures],per_fit["lead_days"],mean),paste0(variant,"_selection_by_lead"))
    fam <- as.data.frame(table(sel$location,sel$variable));names(fam) <- c("location","variable","nonzero")
    fam <- fam[fam$nonzero>0,];fam$mean_per_fit <- fam$nonzero/225
    write(fam,paste0(variant,"_selected_families"))
  }
  for(v in key) selection(v)
  invisible(rows)
}
if(sys.nframe()==0L) summarize_results()
