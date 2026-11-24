.weather_preprocess_dir <- local({
  # source() records the script path. With chdir=TRUE it already entered its folder.
  sources <- Filter(function(f) is.character(f$ofile) && length(f$ofile)==1L, sys.frames())
  if (length(sources)) {
    caller <- sources[[length(sources)]]
    script <- if (is.null(caller$owd)) caller$ofile else basename(caller$ofile)
  } else {
    args <- commandArgs(trailingOnly=FALSE)
    file_arg <- grep("^--file=", args, value=TRUE)
    script <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else args[match("-f", args)+1L]
    if (is.na(script)) script <- "preprocess_weather.R"
    if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed=TRUE)
  }
  dirname(normalizePath(script, mustWork=TRUE))
})

season_basis <- function(date, start, pairs=2L) {
  phase <- 2*pi*as.numeric(as.Date(date)-as.Date(start))/365.2425
  if(pairs==0L) return(matrix(1,length(phase),1L))
  cbind(1, sin(outer(phase,seq_len(pairs))), cos(outer(phase,seq_len(pairs))))
}

learn_scaling <- function(x) {
  med <- apply(x,2,median,na.rm=TRUE)
  med[!is.finite(med)] <- 0
  for (j in seq_len(ncol(x))) x[!is.finite(x[,j]),j] <- med[j]
  sd0 <- apply(x,2,sd)
  sd0[!is.finite(sd0) | sd0==0] <- 1
  list(median=med, center=colMeans(x), sd=sd0)
}

apply_scaling <- function(x, p) {
  for (j in seq_len(ncol(x))) x[!is.finite(x[,j]),j] <- p$median[j]
  sweep(sweep(x,2,p$center,"-"),2,p$sd,"/")
}

season_coefficients <- function(x,basis,train) {
  apply(x[train,,drop=FALSE],2,function(y) {
    ok <- is.finite(y)
    lm.fit(basis[train[ok],,drop=FALSE],y[ok])$coefficients
  })
}

# Source this file to prepare model inputs, or run it to export a fitted model's matrices.
target_fields <- c("Max","Mean","Min")
temperature_cols <- function(a) paste0(a,"_",target_fields,".TemperatureF")
# An origin is the last observed input day. All five target dates must exist
# within the requested period. Observations preceding its first day can be used.
forecast_origins <- function(dates,years) {
  if(length(dates)<10L) return(integer())
  in_period <- as.integer(format(dates,"%Y")) %in% years
  candidates <- seq.int(5L,length(dates)-5L)
  ok <- Reduce(`&`,lapply(1:5,function(h) in_period[candidates+h]))
  candidates[ok]
}

result_directory <- function(code_directory, remote_types="all") {
  variant <- switch(remote_types, mean="mean_only", all="baseline")
  file.path(dirname(code_directory),"main_results",variant)
}

weather_numeric <- function(weather, stations, field, low, high) {
  part <- weather[,paste0(stations,"_",field),drop=FALSE]
  if (!all(vapply(part,is.numeric,logical(1)))) {
    part[] <- lapply(part,function(x) suppressWarnings(as.numeric(as.character(x))))
  }
  z <- as.matrix(part)
  storage.mode(z) <- "double"
  z[!is.finite(z) | z<low | z>high] <- NA_real_
  colnames(z) <- stations
  z
}

clean_weather <- function(weather, stations) {
  dates <- as.Date(rownames(weather))
  stopifnot(length(dates)>=5L, !anyNA(dates),all(diff(dates)==1),
            !anyDuplicated(stations))
  parts <- lapply(target_fields,function(k) weather_numeric(weather,stations,paste0(k,".TemperatureF"),-100,150))
  bad <- which(parts[[1]]<parts[[2]] | parts[[2]]<parts[[3]] | parts[[1]]<parts[[3]])
  cn <- unlist(lapply(stations,temperature_cols),use.names=FALSE)
  temp <- matrix(NA_real_,nrow(weather),length(cn),dimnames=list(rownames(weather),cn))
  for(j in seq_along(parts)) {
    parts[[j]][bad] <- NA_real_
    temp[,seq(j,length(cn),by=3)] <- parts[[j]]
  }
  pressure <- weather_numeric(weather,stations,"Mean.Sea.Level.PressureIn",20,35)
  list(dates=dates,stations=stations,temp=temp,pressure=pressure)
}

load_weather_project <- function(directory) {
  required <- c(file.path(dirname(directory),"data","weather.rda"),file.path(directory,"resources","coordinates.csv"))
  if(!all(file.exists(required))) stop("Missing project input: ",paste(required[!file.exists(required)],collapse=", "),call.=FALSE)
  data <- new.env(parent=emptyenv())
  load(required[1],envir=data)
  coords <- read.csv(required[2])
  stopifnot(nrow(coords)==98L,!anyDuplicated(coords$code))
  targets <- strsplit("ABQ ATL BOS DEN DFW DLH DTW EYW IAD LAN LAX SEA SFO SLC TVC"," ")[[1]]
  list(raw=data$weather,coords=coords,targets=targets,obs=clean_weather(data$weather,coords$code))
}

airport_distances <- function(coords,a) {
  lat <- coords$lat*pi/180
  lon <- coords$lon*pi/180
  i <- match(a,coords$code)
  12742*asin(sqrt(sin((lat-lat[i])/2)^2+cos(lat)*cos(lat[i])*sin((lon-lon[i])/2)^2))
}

five_days <- function(x) {
  blocks <- lapply(4:0,function(lag) {
    index <- seq_len(nrow(x))-lag
    x[replace(index,index<1L,NA_integer_),,drop=FALSE]
  })
  do.call(cbind,blocks)
}

prepare_days <- function(obs,train=NULL,fitted=NULL,pairs=2L) {
  start <- if(is.null(fitted)) obs$dates[1] else fitted$season_start
  basis <- season_basis(obs$dates,start,pairs)
  if(!is.null(fitted)) {
    tc <- fitted$temp_season
    ps <- fitted$pressure_scaling
  } else {
    tc <- if(pairs==0L) matrix(0,1,ncol(obs$temp),dimnames=list(NULL,colnames(obs$temp))) else
      season_coefficients(obs$temp,basis,train)
    ps <- learn_scaling(obs$pressure[train,,drop=FALSE])
  }
  normal <- basis %*% tc
  list(normal=normal,anomaly=obs$temp-normal,local_pressure=obs$pressure,
       pressure=apply_scaling(obs$pressure,ps),temp_season=tc,pressure_scaling=ps,stations=obs$stations)
}

pressure_inputs <- function(p) cbind(p,p-rbind(NA_real_,p[-nrow(p),,drop=FALSE]))

airport_inputs <- function(daily,a,remote_types=c("all","mean")) {
  remote_types <- match.arg(remote_types)
  ext <- setdiff(daily$stations,a)
  local <- cbind(daily$anomaly[,temperature_cols(a)],daily$local_pressure[,a])
  colnames(local) <- paste0(a,c("_Tmax_anomaly","_Tmean_anomaly","_Tmin_anomaly","_Pressure"))
  fields <- if(remote_types=="mean") "Mean" else target_fields
  temperatures <- do.call(cbind,lapply(fields,function(k) {
    x <- daily$anomaly[,paste0(ext,"_",k,".TemperatureF"),drop=FALSE]
    colnames(x) <- paste0(ext,"_T",tolower(k),"_anomaly")
    x
  }))
  p <- daily$pressure[,ext,drop=FALSE]
  raw <- cbind(five_days(local),five_days(temperatures),pressure_inputs(p))
  lags <- function(x) as.vector(outer(colnames(x),paste0("lag",4:0),paste,sep="_"))
  colnames(raw) <- c(lags(local),lags(temperatures),paste0(ext,"_Pressure_lag0"),paste0(ext,"_Pressure_change_lag0"))
  stopifnot(!anyDuplicated(colnames(raw)))
  raw
}

prepare_forecast <- function(model,five_rows,airport) {
  stopifnot(nrow(five_rows)==5L,airport %in% names(model$fits))
  obs <- clean_weather(five_rows,model$stations)
  daily <- prepare_days(obs,fitted=model)
  origin <- obs$dates[5]
  future <- season_basis(origin+1:5,model$season_start) %*% model$temp_season[,temperature_cols(airport),drop=FALSE]
  colnames(future) <- target_fields
  list(raw=airport_inputs(daily,airport,model$remote_types)[5,,drop=FALSE],origin=origin,normal=future)
}

# Decode already-standardized columns. Each column stores a table of completed
# numeric values and integer row indices; this avoids repeating identical data.
get_model_data <- function(data, airport, target, lead, split="train") {
  target <- match.arg(target,c("Min","Mean","Max"))
  split <- match.arg(split,c("train","2013","2014","evaluation","all"))
  key <- paste(airport,target,lead,sep="_")
  case <- data$models[[key]]
  if(is.null(case)) stop("Unknown airport, target or lead: ",key)
  o <- case$observations
  keep <- switch(split,train=o$split=="train",evaluation=o$split!="train",
                 all=rep(TRUE,nrow(o)),o$split==split)
  o <- o[keep,,drop=FALSE]
  x <- vapply(case$feature_ids,function(id) {
    code <- data$daily_codes[[data$code_id[id]]][o$origin_row]
    data$standardized_values[[id]][code]
  },numeric(nrow(o)))
  dimnames(x) <- list(as.character(o$origin),case$feature_names)
  o$origin_row <- NULL
  list(x=x,y=o$target_anomaly_f,observations=o)
}

save_processed_data <- function(directory=.weather_preprocess_dir, remote_types="all") {
  artifacts <- file.path(result_directory(directory,remote_types),"artifacts")
  model <- readRDS(file.path(artifacts,"airport_lasso.rds"))
  stopifnot(identical(model$remote_types,remote_types))
  obs <- load_weather_project(directory)$obs
  temp <- obs$temp
  dates <- obs$dates
  year_id <- as.integer(format(dates,"%Y"))
  train <- which(year_id<=2012)
  daily <- prepare_days(obs,train,model)
  stopifnot(identical(daily$temp_season,model$temp_season),
            identical(daily$pressure_scaling,model$pressure_scaling))
  eval <- forecast_origins(dates,2013:2014)
  eligible <- forecast_origins(dates,unique(year_id[train]))
  origins <- c(eligible,eval)
  daily_codes <- raw_lookup <- list()
  index <- new.env(hash=TRUE, parent=emptyenv())
  value_id <- new.env(hash=TRUE, parent=emptyenv())
  limit <- sum(lengths(model$fits))*length(model$fits[[1]][[1]]$beta)
  standardized <- vector("list",limit)
  code_id <- integer(limit)
  models <- list()
  nvalue <- 0L
  for(a in names(model$fits)) {
    raw <- airport_inputs(daily,a,remote_types)[origins,,drop=FALSE]
    # Local pressure is raw; external pressure is standardized before differencing.
    # Give local and external columns separate names in the shared dictionary.
    source_name <- paste0(c(rep("local/",20),rep("external/",ncol(raw)-20L)),colnames(raw))
    term_id <- integer(ncol(raw))
    for(j in seq_len(ncol(raw))) {
      term <- source_name[j]
      id <- index[[term]]
      if(is.null(id)) {
        id <- length(daily_codes)+1L
        index[[term]] <- id
        v <- unname(raw[,j])
        v[!is.finite(v)] <- NA_real_
        u <- unique(v)
        daily_codes[[id]] <- match(v,u)
        raw_lookup[[id]] <- u
      }
      term_id[j] <- id
    }
    for(f in model$fits[[a]]) {
      sc <- f$scaling
      stopifnot(identical(names(f$beta),colnames(raw)))
      # Hexadecimal floating-point keys identify exact doubles.
      keys <- paste(term_id,sprintf("%a",sc$median),sprintf("%a",sc$center),sprintf("%a",sc$sd),sep=":")
      ids <- integer(length(keys))
      for(j in seq_along(keys)) {
        id <- value_id[[keys[j]]]
        if(is.null(id)) {
          nvalue <- nvalue+1L
          id <- nvalue
          value_id[[keys[j]]] <- id
          v <- raw_lookup[[term_id[j]]]
          v[!is.finite(v)] <- sc$median[j]
          standardized[[id]] <- (v-sc$center[j])/sc$sd[j]
          code_id[id] <- term_id[j]
        }
        ids[j] <- id
      }
      target_col <- paste0(a,"_",f$target,".TemperatureF")
      h <- f$lead
      train_origins <- eligible
      train_origins <- train_origins[is.finite(daily$anomaly[train_origins+h,target_col])]
      stopifnot(length(train_origins)==f$training_n)
      z <- c(train_origins,eval)
      records <- data.frame(origin_row=match(z,origins),origin=dates[z],forecast_date=dates[z+h],
                        split=c(rep("train",length(train_origins)),as.character(year_id[eval+h])),
                        observed_f=unname(temp[z+h,target_col]),seasonal_normal_f=unname(daily$normal[z+h,target_col]),
                        target_anomaly_f=unname(daily$anomaly[z+h,target_col]))
      models[[paste(a,f$target,h,sep="_")]] <- list(feature_ids=ids,feature_names=colnames(raw),observations=records)
    }
  }
  data <- list(training_end=as.Date("2012-12-31"),origins=dates[origins],models=models,daily_codes=daily_codes,
               code_id=code_id[seq_len(nvalue)],standardized_values=standardized[seq_len(nvalue)])
  stopifnot(length(models)==225L)
  saveRDS(data,file.path(artifacts,"weather_preprocessed.rds"),compress=TRUE)
}

if (sys.nframe()==0L) save_processed_data()
