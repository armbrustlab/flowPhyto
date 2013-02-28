summarize <- function(x, channel.clmns=CHANNEL.CLMNS, opp.paths.str='1,2,3') {
  reqd.sds <- c('lat','long','time.UTC','flow.rate', 'bulk.red', 'salinity', 'temperature', 'event.rate', 'fluorescence')
  reqd.clmns <- c('year_day','file','pop', 'evt', 'opp', reqd.sds)
  
  # require(IDPmisc)
    
  reqd.present <- reqd.clmns %in% names(x)
  if(!all(reqd.present))
    stop(paste('you are missing one ore more of the required sds/log or pop columns. Namely: ',paste(reqd.clmns[!reqd.present],collapse=',')))
  file.ct <- length(strsplit(opp.paths.str,',')[[1]])
  
  time.diff <- 3 * file.ct  # 3 minute per file approximation for now!

  if(!'vol.correct' %in% names(x))
    x$vol.correct <- 1
   
  out <- list()
  for(p in unique(x$pop)){
  
    xpp <- subset(x, pop == p)

    meta.df <- 
           data.frame(
	      day = xpp$year_day[1],
	      file   = round(mean(xpp$file, na.rm=TRUE), 2),
	      pop = paste(unique(xpp$pop), collapse=','),
	      resamp = opp.paths.str,
	      time = mean(xpp$time.UTC, na.rm=TRUE),
	      lat    = round(mean(xpp$lat, na.rm=TRUE),5),
	      long   = round(mean(xpp$long, na.rm=TRUE),5),
	      flow = round(mean(xpp$flow.rate, na.rm=TRUE), 3),
	      bulk_red = round(mean(xpp$bulk.red, na.rm=TRUE), 3),
	      salinity = round(mean(xpp$salinity, na.rm=TRUE), 3),
	      temperature = round(mean(xpp$temperature, na.rm=TRUE), 3),
	      event_rate = round(mean(xpp$event.rate, na.rm=TRUE), 3),
	      fluorescence = round(mean(xpp$fluorescence, na.rm=TRUE), 3),
          opp.vol.correction = round(mean(xpp$opp/xpp$evt, na.rm=TRUE),3),
          evt = round(mean(xpp$evt, na.rm=TRUE),2),
          opp = round(mean(xpp$opp, na.rm=TRUE),2),
	      n = nrow(xpp)	    
	   )
    meta.df$conc <- with(meta.df, round(n/(flow * time.diff * opp.vol.correction),4))  # ~Vol/100 correction for OPP filtration
    meta.df$opp.vol.correction <- NULL
    
    for(c in channel.clmns){
            if(nrow(xpp) > 2){
            d <- peaks(density(xpp[,c],from=0, to=2^16, n=2^16))
            meta.df[,paste(c,"_mean",sep="")] <- round(mean(xpp[,c], na.rm=T),3)
            meta.df[,paste(c,"_median",sep="")] <- round(median(xpp[,c], na.rm=T),3)
            meta.df[,paste(c,"_sd",sep="")] <- round(sd(xpp[,c], na.rm=T),3)
            meta.df[,paste(c,"_mode",sep="")] <- round(max(d$x),3)
            meta.df[,paste(c,"_width",sep="")] <- round(d[which(d$x == max(d$x)),"w"],3)
            meta.df[,paste(c,"_npeaks",sep="")] <- nrow(d)  
                }else{
            meta.df[,paste(c,"_mean",sep="")] <- meta.df[,paste(c,"_median",sep="")] <- meta.df[,paste(c,"_sd",sep="")] <- meta.df[,paste(c,"_mode",sep="")] <- meta.df[,paste(c,"_width",sep="")] <- meta.df[,paste(c,"_npeaks",sep="")] <- NA 
               }
            }
 
    # channels <- list()
    # for(stat in c('median','mean', 'sd'))
      # channels[[stat]] <- data.frame(lapply(channel.clmns, function(chnl) 
  				# round(match.fun(stat)(xpp[,chnl], na.rm=TRUE), 1))
  			# )
    # names(channels[['median']]) <- paste(channel.clmns,'_med', sep="")
    # names(channels[['mean']]) <- paste(channel.clmns,'_mean', sep="")
    # names(channels[['sd']]) <- paste(channel.clmns, '_sd', sep='')
     
    # out[[p]] <- cbind.data.frame(meta.df, channels[['median']], channels[['mean']], channels[['sd']])
    
    out[[p]] <- meta.df
    
	}
		
 	do.call(rbind.data.frame, out)
	
}



summarizeFile <- function(opp.paths, pop.names, output.path=getCruisePath(opp.paths[1])){
  
  ## name the opp.paths vector so joinSDS can index it properly
  names(opp.paths) <- sapply(opp.paths, getFileNumber)
  opp.paths.str <- paste(names(opp.paths), collapse=',')
  ##################################################################
  ## CONCATENATE OPP Filtered Evt & Classification Consensus FILES #
  ##################################################################
  filter.df <- do.call(rbind.data.frame, lapply(opp.paths, readSeaflow, add.yearday.file=TRUE))
  consen.df <- do.call(rbind.data.frame, lapply(opp.paths, readConsensusFile)) 
  classed <- cbind.data.frame(filter.df, consen.df)
  
  ## event counts from before and after filtration are needed for proper conc calculations
  nrow.opp <- sapply(opp.paths, function(p) readSeaflow(              p , count.only=TRUE))
  nrow.evt <- sapply(opp.paths, function(p) readSeaflow(sub('.opp','',p), count.only=TRUE))
  #classed$vol.correct <- rep(nrow.opp/nrow.evt, times=nrow.opp)
  classed$opp <- rep(nrow.opp, times=nrow.opp)
  classed$evt <- rep(nrow.evt, times=nrow.opp)
    
  ################################################
  ## WRITE RESAMPLED AGGREGATE STATISTICS FILES ##
  ################################################

  class.pop <- subset(classed, pop %in% pop.names)
  class.jn <- do.call(rbind.data.frame, by(class.pop, list(class.pop$file), joinSDS, opp.paths) )
  stats <- summarize(class.jn, opp.paths.str=opp.paths.str)

  out.paths <- sapply(opp.paths, .createOutputPath, output.path)

  dir.create(dirname(out.paths[length(out.paths)]))
  
  stats.path <- .createStatsPath(out.paths, pop.names)
  write.delim(stats, stats.path)
}


.createStatsPath <- function(paths, pop.names){
  pnames <- sapply(paths, getFileNumber)
  if(length(pnames) == 1)
    path.range = pnames
  else
    path.range <- paste(pnames[1], pnames[length(pnames)],sep=',')  
  paste(dirname(paths[1]),'/', path.range,'.', paste(pop.names,collapse=','),'.stat.tab', sep='')
} 

.combinePopStatFiles <- function(cruise.dir='.'){
  pops <- NULL
  for(this.path in getCruiseFiles(cruise.dir, ext='\\.stat\\.tab')){
      tab <- read.delim(this.path, as.is=TRUE) 
      pops <- rbind.data.frame(pops, tab)
  }
  return(pops)
}

