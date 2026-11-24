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
    if (is.na(script)) script <- "run_comparisons.R"
    if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed=TRUE)
  }
  dirname(normalizePath(script, mustWork=TRUE))
})
source(file.path(.comparison_dir,"airport_lasso.R"),local=TRUE)

# Temperature comparisons change only the representation of external mean
# temperature. Local inputs and all 194 external pressure features stay fixed.
temperature_designs <- c("temp_near5","temp_pc5","temp_pc15")
comparison_specs <- data.frame(
  name=c(rep("season",3),"local","remote_t","nearest10","pc8","individual",temperature_designs),
  pairs=c(0L,2L,3L,rep(2L,8))
)

external_temperature_design <- function(name,daily,airport,coords,train,transform) {
  x <- airport_inputs(daily,airport,"mean")
  ext <- setdiff(daily$stations,airport)
  source <- sub("_.*","",colnames(x))
  remote_temp <- source!=airport & grepl("_Tmean_anomaly_",colnames(x),fixed=TRUE)
  if(name=="temp_near5") {
    km <- airport_distances(coords,airport)[match(ext,coords$code)]
    nearest <- ext[head(order(km,na.last=NA),5L)]
    stopifnot(length(nearest)==5L)
    return(x[,!remote_temp | source %in% nearest,drop=FALSE])
  }
  k <- switch(name,temp_pc5=5L,temp_pc15=15L)
  means <- daily$anomaly[,paste0(ext,"_Mean.TemperatureF"),drop=FALSE]
  if(is.null(transform)) {
    scaling <- learn_scaling(means[train,,drop=FALSE])
    z <- apply_scaling(means[train,,drop=FALSE],scaling)
    rotation <- eigen(crossprod(z),symmetric=TRUE)$vectors[,seq_len(k),drop=FALSE]
    transform <- list(scaling=scaling,rotation=rotation)
  }
  pcs <- apply_scaling(means,transform$scaling) %*% transform$rotation
  lagged <- five_days(pcs)
  colnames(lagged) <- as.vector(outer(paste0("Tmean_PC",seq_len(k)),paste0("lag",4:0),paste,sep="_"))
  design <- cbind(x[,!remote_temp,drop=FALSE],lagged)
  attr(design,"pc_rotation") <- transform
  design
}

comparison_design <- function(name,daily,airport,coords,train=NULL,transform=NULL) {
  if(name %in% temperature_designs)
    return(external_temperature_design(name,daily,airport,coords,train,transform))
  x <- airport_inputs(daily,airport,"all")
  local_temp <- seq_len(20L) %% 4L != 0L
  if(name=="season") return(x[,which(local_temp),drop=FALSE])
  if(name=="local") return(x[,1:20,drop=FALSE])
  if(name=="remote_t") return(x[,1:1475,drop=FALSE])
  if(name=="individual") return(x)
  ext <- setdiff(daily$stations,airport)
  pressure <- daily$pressure[,ext,drop=FALSE]
  if(name=="pc8") {
    directions <- if(is.null(transform)) eigen(crossprod(pressure[train,,drop=FALSE]),symmetric=TRUE)$vectors[,1:8,drop=FALSE] else transform
    p <- pressure %*% directions
    labels <- paste0("PC",1:8)
  } else {
    stopifnot(name=="nearest10")
    km <- airport_distances(coords,airport)[match(ext,coords$code)]
    near <- head(order(km,na.last=NA),10)
    weights <- 1/pmax(km[near],1e-8)
    p <- pressure[,near,drop=FALSE] %*% (weights/sum(weights))
    labels <- "Nearest10"
  }
  out <- pressure_inputs(p)
  colnames(out) <- c(paste0(labels,"_Pressure_lag0"),paste0(labels,"_Pressure_change_lag0"))
  design <- cbind(x[,1:1475],out)
  if(name=="pc8") attr(design,"pc_rotation") <- directions
  design
}

comparison_split <- function(project,end,origins,final=FALSE,chosen=NULL,dest=NULL,spec=comparison_specs) {
  obs <- project$obs
  dates <- obs$dates
  years <- as.integer(format(dates,"%Y"))
  train <- which(years<=end)
  eligible <- forecast_origins(dates,unique(years[train]))
  daily <- lapply(c(0L,2L,3L),function(k) prepare_days(obs,train,pairs=k))
  names(daily) <- c("0","2","3")
  grid <- 2^(3:-5)
  fixed <- 0.0625
  if(final) {
    models <- file.path(dest,"models")
    dir.create(models,showWarnings=FALSE,recursive=TRUE)
    shared <- lapply(daily,function(d) list(temp_season=d$temp_season,pressure_scaling=d$pressure_scaling,season_start=dates[1],stations=obs$stations))
    saveRDS(list(season=shared,coords=project$coords),file.path(models,"preprocessing.rds"))
  }
  airport_run <- function(a) {
    saved <- list()
    result <- list()
    for(s in seq_len(nrow(spec))) {
      name <- spec$name[s]
      day <- daily[[as.character(spec$pairs[s])]]
      raw <- comparison_design(name,day,a,project$coords,train)
      previous_idx <- NULL
      fit_scalers <- list()
      model_fits <- list()
      for(h in 1:5) {
        columns <- temperature_cols(a)
        ytrain <- day$anomaly[eligible+h,columns,drop=FALSE]
        ytest <- day$anomaly[origins+h,columns,drop=FALSE]
        valid <- lapply(1:3,function(j) which(is.finite(ytrain[,j])))
        for(idx in unique(valid)) {
          kinds <- which(vapply(valid,identical,logical(1),idx))
          n <- length(idx)
          rows <- eligible[idx]
          if(!identical(idx,previous_idx)) {
            sc <- learn_scaling(raw[rows,,drop=FALSE])
            x <- apply_scaling(raw[c(rows,origins),,drop=FALSE],sc)
            xt <- x[seq_len(n),,drop=FALSE]
            xv <- x[n+seq_along(origins),,drop=FALSE]
            gram <- crossprod(xt)/n
            if(name=="individual") eig <- eigen(gram,symmetric=TRUE)
            if(final) fit_scalers[[length(fit_scalers)+1L]] <- sc
            previous_idx <- idx
          }
          y <- ytrain[idx,kinds,drop=FALSE]
          mu <- colMeans(y)
          yc <- sweep(y,2,mu)
          p <- ncol(x)
          xy <- crossprod(xt,yc)/n
          if(name=="individual") {
            projected <- crossprod(eig$vectors,xy)
          }
          for(i in seq_along(kinds)) {
            j <- kinds[i]
            lambdas <- if(name=="individual") grid else fixed
            fit <- .C("fit_path",gram,as.double(xy[,i]),as.integer(p),as.double(lambdas),
                      as.integer(length(lambdas)),path=double(p*length(lambdas)),err=double(length(lambdas)))
            stopifnot(max(fit$err)<1.01e-5)
            b <- matrix(fit$path,p)
            labels <- rep(if(name=="season") paste0("local_temp_fourier_",spec$pairs[s]) else name,length(lambdas))
            penalties <- lambdas
            if(name=="individual") {
              ridge <- (eig$vectors %*% (projected[,i]/outer(eig$values,grid,"+")))
              # Reuse the actual fixed-penalty support from the featured estimator.
              scalar <- .C("fit_path",gram,as.double(xy[,i]),as.integer(p),fixed,1L,path=double(p),err=double(1))
              stopifnot(scalar$err<1.01e-5)
              active <- which(scalar$path!=0)
              ols <- numeric(p)
              if(length(active)) ols[active] <- lm.fit(xt[,active,drop=FALSE],yc[,i])$coefficients
              ols[is.na(ols)] <- 0
              b <- cbind(b,ridge,ols)
              labels <- c(labels,rep("ridge",length(grid)),"post_ols")
              penalties <- c(penalties,grid,NA_real_)
            }
            if(final) {
              keep <- if(name=="individual") which((labels=="individual" & penalties==chosen[["individual"]]) |
                         (labels=="ridge" & penalties==chosen[["ridge"]]) | labels=="post_ols") else seq_along(labels)
              for(k in keep) {
                tag <- if(labels[k] %in% c("individual","ridge")) paste0(labels[k],"_tuned") else labels[k]
                key <- paste(tag,target_fields[j],h,sep="/")
                model_fits[[key]] <- list(target=target_fields[j],lead=h,method=tag,
                           beta=setNames(b[,k],colnames(raw)),intercept=unname(mu[i]),scaler=length(fit_scalers),
                           training_n=n,lambda=if(labels[k]=="post_ols") NA_real_ else penalties[k])
              }
            }
            losses <- function(err,split) {
              ok <- is.finite(err[,1])
              data.frame(method=labels,lambda=penalties,airport=a,lead_days=h,target=target_fields[j],
                         fit_end=end,split=split,n=sum(ok),sse=colSums(err[ok,,drop=FALSE]^2))
            }
            e <- sweep(xv %*% b,1,ytest[,j]-mu[i],"-")
            if(final) {
              et <- sweep(xt %*% b,1,yc[,i],"-")
              z <- list(losses(et,"training"))
              for(yr in 2013:2014) z[[length(z)+1L]] <- losses(e[years[origins+h]==yr,,drop=FALSE],as.character(yr))
              z <- do.call(rbind,z)
            } else z <- losses(e,"cv")
            result[[length(result)+1L]] <- z
          }
        }
      }
      if(final) {
        saved[[paste(name,spec$pairs[s],sep="_")]] <- list(design=name,pairs=spec$pairs[s],
                    pc_rotation=attr(raw,"pc_rotation"),scalers=fit_scalers,fits=model_fits)
      }
    }
    if(final) saveRDS(saved,file.path(models,paste0(a,".rds")))
    do.call(rbind,result)
  }
  workers <- if(.Platform$OS.type=="windows") 1L else as.integer(Sys.getenv("AIRPORT_WORKERS","4"))
  out <- parallel::mclapply(project$targets,airport_run,mc.cores=workers,mc.preschedule=FALSE)
  stopifnot(all(vapply(out,is.data.frame,logical(1))))
  do.call(rbind,out)
}

run_comparisons <- function(directory=.comparison_dir) {
  project <- load_weather_project(directory)
  load_solver(directory)
  dest <- file.path(dirname(directory),"main_results","experiments")
  dir.create(dest,showWarnings=FALSE,recursive=TRUE)
  dates <- project$obs$dates
  cv <- forecast_origins(dates,2008:2012)
  issue <- as.integer(format(dates[cv+1L],"%Y"))
  chunks <- lapply(2007:2011,function(end) {
    message("Comparing designs for CV year ",end+1L)
    comparison_split(project,end,cv[issue==end+1L])
  })
  dev <- do.call(rbind,chunks)
  tuned <- aggregate(cbind(n,sse)~method+lambda,subset(dev,method %in% c("individual","ridge")),sum)
  chosen <- vapply(split(tuned,tuned$method),function(x) x$lambda[which.min(x$sse/x$n)],numeric(1))
  message("Final fits using development-selected penalties: ",paste(names(chosen),chosen,collapse=", "))
  chunks[[6]] <- comparison_split(project,2012,forecast_origins(dates,2013:2014),TRUE,chosen,dest)
  details <- do.call(rbind,chunks)
  write.csv(details,file.path(dest,"loss_by_fit_fold_and_penalty.csv"),row.names=FALSE)
  invisible(details)
}
if(sys.nframe()==0L) run_comparisons()

# Call after sourcing this file. Model names correspond to the score-table IDs.
forecast_comparison <- function(five_rows,airport,method,directory=.comparison_dir) {
  path <- file.path(dirname(directory),"main_results","experiments","models")
  shared <- readRDS(file.path(path,"preprocessing.rds"))
  sets <- readRDS(file.path(path,paste0(airport,".rds")))
  for(set in sets) {
    fits <- Filter(function(f) identical(f$method,method),set$fits)
    if(!length(fits)) next
    saved <- shared$season[[as.character(set$pairs)]]
    obs <- clean_weather(five_rows,saved$stations)
    stopifnot(nrow(five_rows)==5L)
    daily <- prepare_days(obs,fitted=saved,pairs=set$pairs)
    raw <- comparison_design(set$design,daily,airport,shared$coords,transform=set$pc_rotation)[5,,drop=FALSE]
    origin <- obs$dates[5]
    normal <- season_basis(origin+1:5,saved$season_start,set$pairs) %*% saved$temp_season[,temperature_cols(airport),drop=FALSE]
    return(do.call(rbind,lapply(fits,function(f) {
      stopifnot(identical(names(f$beta),colnames(raw)))
      x <- apply_scaling(raw,set$scalers[[f$scaler]])
      j <- match(f$target,target_fields)
      data.frame(airport=airport,origin=origin,forecast_date=origin+f$lead,lead_days=f$lead,
                 target=f$target,prediction_f=normal[f$lead,j]+f$intercept+as.vector(x %*% f$beta))
    })))
  }
  stop("Unknown comparison model: ",method)
}
