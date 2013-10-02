census <- function(v, pop.def=POP.DEF){
  x <- table(v)
  ## now match up these names/numbers to the pop.def to make sure we didn't miss anything entirely
  
  pop.names <- c(rownames(pop.def),'x')
  x <- x[match(pop.names,names(x))] 
  names(x) <- pop.names
  x[is.na(x)] <- 0
  x
}

consensus <- function(mtrx, threshold=0.6){
  # returns vector of consensus support counts (named with the actual consensus population)
  
  if(ncol(mtrx) > 1){
  
  tab.list <- apply(mtrx, 1, table)
  
  max.counts <- sapply(tab.list, function(x) x[which.max(x)])
  
  ret <- data.frame(pop=names(max.counts), support=max.counts)
  
  # recode any consensus with counts below the threshold count to NA meaning 'ambiguous'
  levels(ret$pop) <- c(levels(ret$pop), NA)
  ret$pop[ret$support < floor(threshold*ncol(mtrx))] <- NA
  
  return(ret)

  	}else{
  		mtrx$support <- 3
  		ret <- mtrx
 		return(ret)

  		}
 
  }


consensusFile <- function(opp.path, pattern = '.[0-9]+-class.vct$', output.path=paste(.createOutputPath(opp.path, getCruisePath(opp.path)),'.consensus.vct',sep='')){
  d <- dirname(opp.path)
  vects <- paste(d,'/',list.files(d, pattern=paste('^',basename(opp.path), pattern ,sep='')),sep='')
  mat <- do.call(cbind,lapply(vects, read.delim))
  
  consensus <- consensus(mat)
  
  ## write out the consensus file (one per event file)
  write.delim(consensus, output.path)
  
  invisible(consensus)
}

censusFile <- function(opp.path, map.margin=2, output.path=getCruisePath(opp.path), def.path=paste(getCruisePath(opp.path),'pop.def.tab',sep='')){
  ## opp.path <- paste(REPO.PATH,'/Thompson_0/2009_311/3.evt.opp',sep='')
 
  ######################################################
  ## CONCATENATE Classification FILES into a CONSENSUS #
  ######################################################
  consen.df <- consensusFile(opp.path) 
  
  out.path <- .createOutputPath(opp.path, output.path)
  
  if(!is.na(def.path))
    pop.def <- readPopDef(def.path)
  else
    pop.def <- POP.DEF
 
  ## perform a cross tabulation of the consensus results
  census <- census(consen.df$pop, pop.def=pop.def)
    
  ## and write out the counts to a file (one line per file)
  census.path <- paste(opp.path,'.census.tab', sep='')
  col.ct <- length(census)+2
  write(c('yearday','file',names(census)), file=census.path, ncolumns=col.ct, sep='\t')
  write(c(.getYearDay(opp.path), getFileNumber(opp.path), census), file=census.path, ncolumns=col.ct, sep='\t', append=TRUE)
 

  ##########
  ## PLOT ##
  ##########
  opp <- readSeaflow(opp.path)
  opp$pop <- consen.df$pop
  

  bitmap(paste(out.path, '.class.gif', sep=''), width=.plot.width, height=3/2*.plot.width)

	hist1 <- hist(opp$fsc_small, breaks=seq(0,2^16, by=2^16/25), plot=FALSE)
	hist2 <- hist(opp$chl_small, breaks=seq(0,2^16, by=2^16/25), plot=FALSE)
	hist3 <- hist(opp$pe, breaks=seq(0,2^16, by=2^16/25), plot=FALSE)
	hist4 <- hist(opp$chl_big, breaks=seq(0,2^16, by=2^16/25), plot=FALSE)
	hist5 <- hist(opp$fsc_perp, breaks=seq(0,2^16, by=2^16/25), plot=FALSE)

	def.par <- par(no.readonly = TRUE) # save default, for resetting...
	nf <- layout(matrix(c(2,0,5,0,1,3,4,6,8,0,11,0,7,9,10,12,14,0,16,16,13,15,16,16),6,4,byrow=TRUE), c(3,1,3,1,3), c(1,3,1,3,1,3), TRUE)
	
	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'chl_small', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist2$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'pe', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist3$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'chl_big', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist4$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'chl_small', 'pe', pop.def=pop.def)
	par(mar=c(0,6,1,1))
	barplot(hist2$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist3$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

	par(mar=c(6,6,1,1))
	plotCytogram(opp, 'fsc_small', 'fsc_perp', pop.def=pop.def, add.legend=TRUE)
	par(mar=c(0,6,1,1))
	barplot(hist1$counts, axes=FALSE, space=0, col=NA)
	par(mar=c(6,0,1,1))
	barplot(hist5$counts, axes=FALSE, space=0, horiz=TRUE, col=NA)

  
    ## now add a map plot so we can know where we are
    latlong <- .getSDSlatlon(.getSDS(opp.path))
    if(!is.null(latlong) & !any(sapply(latlong, is.na))){
      cruise.track <- try(.makeSDSlatlonDF(read.delim(paste(getCruisePath(opp.path),'/sds.tab',sep=''))), silent=T) 
      par(mar=c(6,6,1,1))
      plotLatLongMap(latlong[1], latlong[2], cruise.track, margin=as.numeric(map.margin), pch=3, cex=2, lwd=2, col='red', main=paste("file: ",.getYearDay(opp.path),'/',basename(opp.path), sep=""))
      }
   dev.off()
   
par(def.par)   	  

}


###########################
##  Resampling Functions ##
###########################


readConsensusFile <- function(path){
   read.delim(paste(path,'.consensus.vct',sep=''))
}



.sampleFileCounts <- function(x, resample.min=500, resamp.concat.max=5, ydf.names){
   ## group a populations files by their collective ability to sum up to 500
   fi.fact.idx <- 1
   fi.fact.evt.sum <- 0
   fi.fact <- rep(NA, length(x))
   for(i in seq(along=x)){
       fi.fact[i] <- fi.fact.idx
       fi.fact.evt.sum <- fi.fact.evt.sum + x[i]
       if(fi.fact.evt.sum >= resample.min){
	  fi.fact.evt.sum <- 0
   	  fi.fact.idx <- fi.fact.idx + 1
       }
   }
   
   ## filter (via NA) those resample attempts that would require too many concatenations
   fi.fact.tab <- table(fi.fact)
   fi.fact[as.character(fi.fact) %in% names(fi.fact.tab)[fi.fact.tab > resamp.concat.max]] <- NA
   names(fi.fact) <- ydf.names
   return(fi.fact)
}

.resampleRefactorReduction <- function(rsl){
  
  rfctr <- lapply(rsl, function(r) by(names(r), r, function(y) paste(y, collapse=',')))
  reduct <- list()
  for(pop in names(rfctr))
    reduct[[pop]] <- .nv(rep(pop,length(rfctr[[pop]])), rfctr[[pop]])
  flat <- unlist(reduct)
  resamp.fact <- sapply(strsplit(names(flat), '\\.'), '[', 2)
  return(
         do.call(c, list(by(flat, resamp.fact, function(y) paste(y, collapse=','))))
         )
}

createResamplingScheme <- function(cruise.path, resample.min=500, resamp.concat.max = 5){
  census.path <- paste(cruise.path,'/census.tab',sep='')
  census <- read.delim(census.path, as.is=TRUE)
  census <- census[order(census$yearday, census$file),]

  #census <- sapply(census , function(x) x[x< 3] <- 0 } # if there are fewer than 3 cells don't resample
  ydf.names <- paste(census$yearday,'/',census$file,sep='')
 

  o <- lapply(census[,-match(c('yearday','file'),names(census))], .sampleFileCounts, resample.min, resamp.concat.max, ydf.names)
  o <- o[names(o)!='x'] #don't calculate/resample the undefined events
  return(
         .resampleRefactorReduction(o)
         )  
}

combineCensusFiles <- function(cruise.dir='.'){
  census <- NULL
  for(this.path in getCruiseFiles(cruise.dir, ext='\\.census\\.tab')){
      tab <- read.delim(this.path, as.is=TRUE) 
      census <- rbind.data.frame(census, tab)
  }
  return(census)
}


