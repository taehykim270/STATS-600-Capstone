# Input is the loss totals for each training fit and the individual CV and test forecasts.
# You can rescore a saved run without fitting any regressions.
training_loss <- function(observed,predicted) {
  stopifnot(length(observed)>0,length(observed)==length(predicted),all(is.finite(observed)),all(is.finite(predicted)))
  data.frame(n=length(observed),sse=sum((observed-predicted)^2))
}

score_mse <- function(run) {
  training <- run$training
  training$split <- "training"
  p <- rbind(transform(run$cv,split="cv"),transform(run$evaluation,split=format(origin+lead_days,"%Y")))
  stopifnot(!anyDuplicated(p[c("split","origin","airport","lead_days","target")]),
            all(is.finite(p$prediction_f)))
  # Keep every issued forecast in the saved run; only observed targets have a loss.
  p <- p[is.finite(p$observed_f),,drop=FALSE]
  p$sse <- (p$observed_f-p$prediction_f)^2
  p$n <- 1L
  cols <- c("lead_days","target","split","sse","n")
  all <- rbind(training[cols],p[cols],transform(p[p$split %in% c("2013","2014"),cols],split="pooled"))
  by_group <- function(keys) {
    sums <- aggregate(all[c("sse","n")],all[c(keys,"split")],sum)
    if(identical(keys,"target")) {
      total <- aggregate(sums[c("sse","n")],list(split=sums$split),sum)
      total$target <- "All"
      sums <- rbind(sums,total)
    }
    out <- unique(sums[keys])
    for(s in c("training","cv","2013","2014","pooled")) {
      z <- sums[sums$split==s,]
      label <- if(s %in% c("2013","2014")) paste0("mse_",s) else s
      z[[label]] <- z$sse/z$n
      z[[paste0(label,"_n")]] <- z$n
      out <- merge(out,z[c(keys,label,paste0(label,"_n"))],by=keys)
    }
    out$pooled_rmse_f <- sqrt(out$pooled)
    order_target <- match(out$target,c("Min","Mean","Max","All"))
    if("lead_days" %in% keys) {
      out <- out[order(out$lead_days,order_target),]
    } else out <- out[order(order_target),]
    rownames(out) <- NULL
    out
  }
  test <- p[p$split %in% c("2013","2014"),]
  test$forecast_date <- test$origin+test$lead_days
  calendar <- aggregate(cbind(sse,n)~forecast_date+lead_days+target,test,sum)
  calendar$mse <- calendar$sse/calendar$n
  list(summary=by_group("target"),by_lead=by_group(c("lead_days","target")),
       by_date=calendar[order(calendar$forecast_date,calendar$target),c("forecast_date","lead_days","target","n","mse")])
}

score_saved_predictions <- function(path) {
  run <- readRDS(path)
  s <- score_mse(run)
  output <- if(basename(dirname(path))=="artifacts") dirname(dirname(path)) else dirname(path)
  details <- file.path(output,"details")
  dir.create(details,recursive=TRUE,showWarnings=FALSE)
  write.csv(s$summary,file.path(output,"temperature_scores.csv"),row.names=FALSE)
  write.csv(s$by_lead,file.path(output,"mse_by_lead_target.csv"),row.names=FALSE)
  write.csv(s$by_date,file.path(details,"mse_by_calendar_date.csv"),row.names=FALSE)
  issued <- run$evaluation
  issued$n_scored <- as.integer(is.finite(issued$observed_f))
  issued$sse <- ifelse(issued$n_scored==1L,(issued$prediction_f-issued$observed_f)^2,0)
  issued$n_issued <- 1L
  windows <- aggregate(cbind(sse,n_scored,n_issued)~origin+airport+target,issued,sum)
  stopifnot(all(windows$n_issued==5L))
  windows$mse <- ifelse(windows$n_scored>0L,windows$sse/windows$n_scored,NA_real_)
  windows <- windows[order(windows$origin,windows$airport,windows$target),]
  write.csv(windows,file.path(details,"mse_by_window.csv"),row.names=FALSE)
  s$by_window <- windows
  s
}

if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=1L) stop("Usage: Rscript score_mse.R path/to/forecast_errors_input.rds",call.=FALSE)
  print(score_saved_predictions(args[1])$by_lead)
}
