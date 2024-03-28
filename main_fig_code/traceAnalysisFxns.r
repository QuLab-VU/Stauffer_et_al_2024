require(tidyverse)

plotTraceMean <- function(ids, dat, autorange_y=FALSE, yl=c(500,2000)) {
    dat <- dat[dat$reg_id %in% ids,]
    
    min_y <- min(dat$intensity_mean)
    max_y <- max(dat$intensity_mean)
    
    if(autorange_y)  yl <- c(min_y,max_y) 

    p <- ggplot(dat) +
        geom_line(aes(x = im_idx, y = intensity_mean)) +
        facet_grid(rows = vars(reg_id)) + 
        theme(strip.text.y = element_text(angle = 0)) + 
        ylim(yl)
    p
}

plotTraceFura <- function(ids, dat, yl=c(0,1) ) plotTraceMean(ids=ids, dat=dat, yl=yl)

plotTraceNorm <- function(ids, dat, yl=c(0,2500), tit=NA) {
    dat <- dat[dat$reg_id %in% ids,]
    p <- ggplot(dat) +
        geom_line(aes(x = im_idx, y = int_mean_norm)) +
        ylim(yl) +
        facet_grid(rows = vars(reg_id)) + 
        theme(strip.text.y = element_text(angle = 0)) + 
        ggtitle(tit)
    p
}

plotTraceWellNorm <- function(ids, dat, yl=c(0,2500), tit=NA) {
    dat <- dat[dat$reg_id %in% ids,]
    p <- ggplot(dat) +
        geom_line(aes(x = im_idx, y = int_well_mean_norm)) +
        ylim(yl) +
        facet_grid(rows = vars(reg_id)) + 
        theme(strip.text.y = element_text(angle = 0)) + 
        ggtitle(tit)
    p
}

peakFinder <- function(y,lag,threshold,influence) {
  # Empirically determined optimal parameters: lag=25, threshold=4, influence=0.01  
  signals <- rep(0,length(y))
  filteredY <- y[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  avgFilter[lag] <- mean(y[0:lag], na.rm=TRUE)
  stdFilter[lag] <- sd(y[0:lag], na.rm=TRUE)
  for (i in (lag+1):length(y)){
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
      if (y[i] > avgFilter[i-1]) {
        signals[i] <- 1;
      } else {
        signals[i] <- -1;
      }
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0
      filteredY[i] <- y[i]
    }
    avgFilter[i] <- mean(filteredY[(i-lag):i], na.rm=TRUE)
    stdFilter[i] <- sd(filteredY[(i-lag):i], na.rm=TRUE)
  }
  return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter))
}

getMeanInt <- function(vals, coln, idx) {
    a <- matrix(vals,nrow=length(idx))
    rownames(a) <- idx
    colnames(a) <- coln
    mean_signal <- apply(a,1,mean)
    
    norm_vals <- as.vector(sapply(coln, function(cn) {
        out <- a[,cn] - mean_signal
        out - min(out)
        }))
    norm_vals    
}

prepIntValData <- function(filepath) {
    d <- read.csv(filepath)
    # drop first column (merged rownames from Python); use label to identify each object over time
    d <- d[,-1]
    # rename colname label to reg_id to work with plotTrace functions
    colnames(d)[colnames(d)=="label"] <- "reg_id" 
    
    obj_info <- d[!duplicated(d$reg_id),c("reg_id","centroid.0","centroid.1","area")]
    d <- d[,c("reg_id","intensity_mean")]
    
    reg_ids <- unique(d$reg_id)
    keep_ids <- obj_info[obj_info$area > 100,"reg_id"]
    
    n_obj <- nrow(obj_info)
    
    d$im_idx <- rep(1:(nrow(d)/n_obj), each=n_obj)
    
    out <- list(well=wellFromPath(filepath),
                reg_info=obj_info,
                intvals=d,
                keep_ids=keep_ids)
    
    return(out)
}

wellFromPath <- function(filepaths) gsub("Well","",sapply(strsplit(filepaths, "_"), "[[", 2))


