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
    if (is.na(script)) script <- "plot_selected_airports.R"
    if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed=TRUE)
  }
  dirname(normalizePath(script, mustWork=TRUE))
})
source(file.path(.comparison_dir,"airport_lasso.R"),local=TRUE)

plot_selected_airports <- function(directory=.comparison_dir) {
  root <- file.path(dirname(directory),"main_results")
  dest <- file.path(root,"baseline","maps")
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  boundaries <- read.csv(file.path(directory,"resources","state_outlines.csv.gz"))
  places <- read.csv(file.path(directory,"resources","airport_places.csv"),na.strings=c("","NA"))
  coord <- read.csv(file.path(directory,"resources","coordinates.csv"))
  sel <- read.csv(file.path(root,"baseline","artifacts","selected_predictors.csv"))
  sel <- subset(sel,location=="external" & target=="Mean")
  project <- function(lon,lat) {
    r <- pi/180
    n <- (sin(29.5*r)+sin(45.5*r))/2
    rho <- sqrt(cos(29.5*r)^2+2*n*(sin(29.5*r)-sin(lat*r)))/n
    angle <- n*(lon+96)*r
    cbind(x=rho*sin(angle),y=-rho*cos(angle))
  }
  bg <- project(boundaries$lon,boundaries$lat)
  parts <- split(seq_len(nrow(bg)),interaction(boundaries$state,boundaries$part,drop=TRUE))
  xy <- project(coord$lon,coord$lat)
  rownames(xy) <- coord$code
  for(a in c("DTW","ATL","SFO")) {
    city <- places[match(a,places$code),]
    png(file.path(dest,paste0("locations_",a,".png")),width=1650,height=1120,res=145,bg="white")
    par(mfrow=c(2,3),mar=c(.9,.2,2,.2),oma=c(2,1,4,1))
    for(h in 1:5) {
      z <- subset(sel,target_airport==a & lead_days==h)
      groups <- split(z,z$source_airport)
      col <- vapply(groups,function(v) {
        t <- any(grepl("^T",v$variable))
        p <- any(grepl("Pressure",v$variable))
        if(t && p) "#794E95" else if(t) "#BB5023" else "#217396"
      },character(1))
      mapped <- names(col)[is.finite(xy[names(col),1])]
      plot(bg,type="n",asp=1,axes=FALSE,xlab="",ylab="")
      for(part in parts) polygon(bg[part,],col="#F1F3EF",border="#BAC5BE",lwd=.40)
      points(xy[mapped,,drop=FALSE],pch=21,bg=col[mapped],col="white",cex=.9,lwd=.5)
      points(xy[a,1],xy[a,2],pch=23,bg="#163E4D",col="white",cex=1.4)
      title(sprintf("Lead %d   |   %d selected airports",h,length(groups)),cex.main=.88)
      top <- z[which.max(abs(z$coefficient_per_training_sd)),]
      state <- places[match(top$source_airport,places$code),]
      day <- if(top$variable=="Pressure_change") "t minus (t - 1)" else paste0("lag ",top$lag_days)
      mtext(sprintf("Largest |coefficient|: %s (%s), %s",top$source_airport,state$state,day),side=1,line=-.25,cex=.63)
    }
    plot.new()
    legend("center",legend=c("Temperature","Pressure","Both","Forecast airport"),
           pt.bg=c("#BB5023","#217396","#794E95","#163E4D"),pch=c(21,21,21,23),bty="n",cex=.9,title="Selected observations")
    mtext(paste0(city$city,", ",city$state_name," (",a,")"),outer=TRUE,side=3,line=1.8,cex=1.5,font=2)
    mtext("External inputs retained by the five mean-temperature regressions",outer=TRUE,side=3,line=.4,cex=.93)
    mtext("Fitted through 2012. Markers show stations with known coordinates. Boundaries: U.S. Census Bureau (2023). Locations: OpenFlights.",outer=TRUE,side=1,line=.45,cex=.76)
    dev.off()
  }
}
if(sys.nframe()==0L) plot_selected_airports()
