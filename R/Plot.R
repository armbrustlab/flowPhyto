plotStatMap <- function(df, pop, z.param, margin=0.1, zlab=z.param, ma=1, xlim=NULL, ylim=NULL, main=paste(pop,"-population"), track=NULL, cex=1,...){   
	if(missing(pop)){
	  warning("no value found for pop param. using entire stats dataframe")
	  df.pop <- df
	}else{
	  df.pop <- df[df$pop==pop,] # subset the pop of interest
	}
	z <- df.pop[,z.param] 
	z <- z[is.finite(z)]
	try(z <- rollmean(z, k = ma), silent=TRUE)
	w <- na.exclude(z)
		if(is.na(ma)){
		z <- smooth.spline(z, df = length(unique(z)))
		z <- z$y
			}
	plotLatLongMap(long=df.pop[,"long"], lat=df.pop[,"lat"], margin=margin, col=.makeMapColorGradient(z), legend=pretty(quantile(w, ,seq(0,1,by=.25))), zlab=zlab, xlim=xlim, ylim=ylim, main=main,track=track, cex=cex, ...) 
	
	}

.makeMapColorGradient <- function(y, n.breaks = 20, add.legend=TRUE, log.base = NULL){
  y <- y[is.finite(y)]
  breaks <- seq(min(y, na.rm=TRUE), max(y, na.rm=TRUE), length.out = n.breaks)
  #breaks <- pretty(breaks, n=length(breaks)) # commented out because n only give 'desired' n
  if(!is.null(log.base))
    breaks <- log.base^breaks
  break.labels <- as.character(round(breaks,3))
  quantile <- cut(y, n.breaks, labels = break.labels)

  legend.vect <- .rainbow.cols(n.breaks)
  names(legend.vect) <- break.labels

  ## pcts <- as.numeric(quantile)/n.breaks
  ## ret <- gray(pcts)
  ## ret <- hsv(h=rgb2hsv(col2rgb(.pop.colors[this.pop]))['h',], s=pcts, v=1)
  
  ret <- .rainbow.cols(n.breaks)[quantile]
  if(add.legend)
    attributes(ret) <- list(legend=legend.vect)
  return(ret)
}

plotLatLongMap <- function(lat, long, track=NULL, margin=2, col='red', legend=NULL, pch=20, cex=1.5, lwd=1, lty=2, xlim=NULL, ylim=NULL, xlab="Longitude (deg W)",ylab="Latitude (deg N)",zlab=NA, ...){
  ## plot longitude and latitude on a map
  require(maps, quietly=T)
  require(mapdata, quietly=T)
  require(plotrix, quietly=T)
  
  map.type <- 'worldHires'
  
    if(is.null(track) | class(track)=='try-error')
 	track <- data.frame(lat=lat, long=long)
    
    if(is.null(xlim))
    	    xlim <- c(min(track$long, na.rm=TRUE)-margin, max(track$long, na.rm=TRUE)+margin)    	    
    if(is.null(ylim))
    	    ylim <- c(min(track$lat , na.rm=TRUE)-margin, max(track$lat , na.rm=TRUE)+margin)
  	
  	if(xlim[1] < 0 & xlim[2] > 0){
  		neg.long <- subset(track, long < 0)
  		track[row.names(neg.long), "long"] <- neg.long$long + 360
  		xlim <- c(min(track$long, na.rm=TRUE)-margin, max(track$long, na.rm=TRUE)+margin)
  		long <- na.exclude(long); lat <- na.exclude(lat)
  		long[long < 0] <- long[long < 0] + 360
  		map.type <- 'world2Hires'
  			}
  
  # plot the cruise track as gray line back-ground
  plot(track$long, track$lat, xlim=xlim, ylim=ylim, lwd=lwd,lty=lty,
       xlab=xlab,ylab=ylab,
       pch=20, col='gray', type='o', asp=1, ...)
  try(maps::map(map.type, fill=FALSE, col='black',add=TRUE))
  points(long, lat, xlim=xlim, ylim=ylim, pch=pch, cex=cex, col=col)
	
  if(!is.null(legend)){
    ylim <- par('usr')[c(3,4)]
    xlim <- par('usr')[c(1,2)]

    color.legend(xlim[2], ylim[1], xlim[2] + diff(xlim)/40, ylim[2], 
    	legend=legend, rect.col=.rainbow.cols(100), gradient='y',align='rb',cex=cex,...)
	mtext(zlab, side=4, line=3,cex=cex)  

  }
  if(length(long) == 1)
    mtext(paste("Long/Lat: ",round(mean(long),3),'/', round(max(lat),3)),col='red', line=-2, )
}


.plotStatsByPop <- function(df, x.var, y.var, labs, colrs=POP.DEF$color,  ...){
  plot( df[,x.var], df[,y.var], pch='', xlab=x.var, ylab=labs[y.var], ...)
  by(df, list(df$pop), function(z) points(z[,x.var], z[,y.var], col=colrs[as.character(z$pop[1])]))
  try(
    by(df, list(df$pop), function(z) lines(smooth.spline(z[,x.var], z[,y.var]), col=colrs[as.character(z$pop[1])]))
    )
  legend(min(df[,x.var]), max(df[,y.var]), col=colrs, legend=names(colrs), pch=1)
}



############################################
## SECTION FOR BATCH RUNNING OF THIS FILE ##
############################################
.plot.args <- c('cruise','x.var','y.vars','date.range','pops') #'hour.range','lat.range','long.range',
.indepnd.vars <- c('map','lat','long','time')
.outcome.vars <- c('conc','fsc','chl')
.y.var.labs <- c("Concentration, 10^6 cells/L","avg FSC [Size]","avg Chlorophyl, rfu")
names(.y.var.labs) <- .outcome.vars
.additional.vars <- c('temperature','salinity')
.plot.width <- 9

plotCruiseStats <- function(cruise, x.var='map', y.vars=c("conc", "chl_small","fsc_small"), pops=c('ultra','synecho'), sds.var=NULL,
           ## hour.range=c(00,24), lat.range=c(-90,90), long.range=c(-180,180),  #none of these are used yet
			date.range=as.POSIXct(c('2009-01-01','2099-12-31'), tz='UTC'), output.path=paste(REPO.PATH, cruise,'/plots/',sep=''), ...){
			
          cruise.dir <- dirname(cruise)
    	  if(cruise.dir == '.' & substr(cruise,1,1)!= '.')	
    	    repo.dir <- paste(REPO.PATH,'/',sep='')
    	  else
    	    repo.dir <- paste(cruise.dir,'/',sep='')
    	  
    	  cruise <- basename(cruise)
    	
	  query.db <- FALSE
	  if(query.db){#if query.db is passed to this file then query the database for the file instead of loading it all
	    ranges <- list()
	    ranges['date'] <- date.range
	    pop.stats <- .queryStats(cruise, x.var, y.vars, ranges, pops)
	  }else{
	    ## Read the tab delimited summary file into memory
	    pop.stats <- read.delim(paste(repo.dir, cruise,'/stats.tab',sep=''))
	    pop.stats <- subset(pop.stats, !is.na(resamp))  # because we're now including resamples less than 500
	    ## convert the character representation of time into R's POSIX format
	    pop.stats$time <- as.POSIXct(pop.stats$time)
	    ## subset it acording to the date and pop parameters
	     pop.stats <- subset(pop.stats, date.range[1] <= time & time <= (date.range[2] + 60*60*24))
	     pop.stats <- subset(pop.stats, pop.stats$pop %in%  pops)

	  }

	     ## figure out the populations used
	     pop.defs <- read.delim(paste(repo.dir, cruise,'/pop.def.tab',sep=''))
	     pops.avail.v <- as.character(pop.defs$abrev)
	     .pop.colors <- as.character(pop.defs$color[pop.defs$abrev %in% pops])  #pop.stats$pop
	     names(.pop.colors) <- pops

	  ## set up the directory and file for plotting
		plot.arg.lst = list()
		for(arg in .plot.args)
		    plot.arg.lst[[arg]] <- paste(get(arg), collapse=',')
		plot.name.string <- paste(plot.arg.lst, collapse='.')

	  plot.dir <- paste(output.path,'/',sep='') # safety check
	  if(!file.exists(plot.dir))
	    dir.create(plot.dir)

	  ## prepare (optional) SDS information
	  plot.add.var <- FALSE
	  if(!is.null(sds.var))
	    if(sds.var %in% .additional.vars)
	      plot.add.var <- TRUE


	    sds <- read.delim(paste(repo.dir, cruise,'/sds.tab',sep=''))
	    ## convert the character representation of time into R's POSIX format
	    sds$time <- as.POSIXct(sds$time)
	    sds$salinity <- sds$SALINITY
	    sds$temperature <- sds$OCEAN.TEMP
	    #lat.lon <- t(sapply(1:nrow(sds), function(x) .getSDSlatlon(sds[x,])))
	    #sds$lat <- lat.lon[,1]
	    #sds$long <- lat.lon[,2]
   	 cruise.track <- try(.makeSDSlatlonDF(read.delim(paste(repo.dir, cruise,'/sds.tab',sep='')))) 

	if(plot.add.var){
	    plot.name.string <- paste(plot.name.string, '.', sds.var, sep='')
	 }


	  ## set up the actual bitmap / gif
	  plot.height <- .75 * .plot.width
	  if(x.var == 'map')
	    plot.height <- plot.height * length(y.vars)
	  bitmap(paste(plot.dir, plot.name.string, '.gif', sep=''), height=plot.height, width=.plot.width)

	  ## generate the actual plot  
	  if(x.var %in% c('lat','long','time')){
	    par(mfrow=c(length(y.vars),1))
	    for(y.var in y.vars){
	      par(mar=c(5, 4, 4, 4) + 0.1)
	      .plotStatsByPop(df=pop.stats, x.var, y.var, labs= .y.var.labs, colrs=.pop.colors, main=paste(cruise,'Seaflow'), ...)
	      if(plot.add.var){
		par(new=TRUE)
		plot(sds[,x.var], sds[,sds.var], type='l', col='gray', xlab='', ylab='', axes=FALSE) #bty="n"
		axis(4) 
		mtext(sds.var, side=4, line=3, cex.lab=1, col='gray80')
	      }

	    }
	  }else if(x.var == 'map'){
	    margin <- .1
	    
	    par(mfrow=c(length(y.vars),length(pops)))
	    single.pop <-  length(pops) == 1 # & length(y.vars) == 1 
	    for(y.var in y.vars){
	      log.yvar <- log(pop.stats[,y.var], 10)

	      for(this.pop in pops){
		is.pop <- as.character(pop.stats$pop)==this.pop

		long <- pop.stats$long
		lat <- pop.stats$lat

		## filter gradient & create a legend only if the plot is big enough to fit it in
		if( single.pop){
		  #xlims <- ylims <- NULL #auto plot area is OK
		  log.yvar <- log.yvar[is.pop]
		}
		col.vect <- .makeMapColorGradient(log.yvar, add.legend=single.pop)
		if(!single.pop){
		  # use the global max and min of lat and long to determine uniform plot areas
		  #xlims <- c(min(long, na.rm=TRUE)-margin, max(long, na.rm=TRUE)+margin)
		 # ylims <- c(min(lat, na.rm=TRUE)-margin, max(lat, na.rm=TRUE)+margin)
		  col.vect <- col.vect[is.pop]
		}
		legend.vect <- attr(col.vect,'legend')  # null for non single plots

		plotLatLongMap(lat=lat[is.pop], long=long[is.pop], track=cruise.track, 
				col=col.vect, legend=legend.vect, margin=.2,
			       main=paste(cruise, 'SeaFlow\n',this.pop), zlab=.y.var.labs[y.var], ...) #, xlims=xlims, ylims=ylims)

	      }
	    }
	  }else{
	    stop("this shouldn't be happening... plot type is funky")
	  }

	  ## close the plotting device
	  dev.off()
}


.prePlotLevel2 <- function(cruise.path, output.path='.', log.dir='.', debug=FALSE){

  n <- 0
  ##.plot.args <- c('cruise','date.range', 'x.var','y.vars','populations')  
  ##.indepnd.vars <- c('map','lat','long','time')
  ##.outcome.vars <- c('conc','fsc','chl')
  plotRpath <- paste(.SOURCE.DIR,'/Plot.R', sep='')
  if(!file.exists(log.dir))
    dir.create(log.dir)

  sds <- read.delim(paste(cruise.path, '/sds.tab', sep=''))
  stats <- read.delim(paste(cruise.path, '/pop.stats.tab', sep=''))

  cruise.name <- basename(cruise.path)
  cmd.0 <- paste("R CMD BATCH --no-save --no-restore --cruise=",cruise.name, sep='')
  cmd.0 <- paste(cmd.0,' --date.range=', paste(range(as.Date(sds$time)), collapse=','), sep='')
  cmd.0 <- paste(cmd.0,' --output.path=', output.path,sep='')

  pops <- as.character(unique(stats$pop))
  # add all possible pop combos
  pops <- c(pops, paste(pops, collapse=','))
  # add popular pop combos (from 2 to 5 in count)
  popular <- table(stats$pop, is.na(stats$resamp), useNA='ifany')[,'FALSE']
  popular <- popular[rev(order(popular))]
  pops <- c(pops, sapply(2:5, function(i) paste(names(popular[1:i]),collapse=',')))
  # add all possible combinations of "important" pops (default for the web interface)
  important <- c('pico','ultra','nano','synecho', 'crypto','cocco','prochloro')
  pops <- c(pops, do.call(c,sapply(2:length(important)-1, function(j) apply(combn(important,j), 2, paste, collapse=','))))
  pops <- c(pops, paste(important, collapse=','))

  .outcome.vars <- c(.outcome.vars, paste(.outcome.vars, collapse=','))
  .outcome.vars <- c(.outcome.vars, apply(combn(.outcome.vars,2), 2, paste, collapse=','))

  add.vars <- c('', .additional.vars)


  
  for(x in .indepnd.vars){
    cmd.x <- paste(cmd.0,' --x.var=', x, sep='')
    for(y in .outcome.vars){      
      cmd.y <- paste(cmd.x,' --y.vars=', y, sep='')
      for(pop in pops){
        cmd.p <- paste(cmd.y,' --populations=', pop, sep='')
        for(add in add.vars){
          if(add != '' & x != 'map')
            cmd.p <- paste(cmd.p, ' --sds.var=',add, sep='')
          cmd <- paste(cmd.p, ' ',plotRpath)          
          n <- n+1
          log.path <- paste(log.path, basename(cruise.path),'.',x,'.',y,'.',pop,'.',add,'.log',sep='')
          if(debug)
            print(cmd)
          else
            system(paste("echo '",cmd, log.dir,"' | qsub"))
        }
      }
    }    
  }
 
  #--x.var=map --y.vars=conc,chl,fsc --date.range=2009-11-07,2009-11-09  --populations=ultra --custom.out.dir=. Plot.R
  return(n)  
}

