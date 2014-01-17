##### NEW FUNCTION !!! ####

filter <- function(events, width=1, notch=1, origin=NA, do.plot=FALSE){
  notch <- as.numeric(notch)
  width<- as.numeric(width)
  origin <- as.numeric(origin)
   
   if(any(max(events[,-c(1,2)]) > 10^3.5)){
    stop(paste("ERROR, data are not LOG-transform"))
   }
   
   ##### FILTRING OPP #####
  		detected <- subset(events, D1 > 1 & D2 > 1) # filtering particles not detected by D1 or D2
  		unsaturated <- subset(detected, D1 < max(events[,"D1"]) & D2 < max(events[,"D1"])) # filtering particles with saturated signals on D1 or D2
		unsaturated[,-c(1,2)] <- (log10(unsaturated[,-c(1,2)])/3.5)*2^16 ## linearize the LOG transformed data

		if(is.na(origin)){
		  origin.unsaturated <- subset(unsaturated, D1 > 10000 & D2 > 10000) # Exclude potenital electrical noise from calculation.	
		  origin <- median(origin.unsaturated$D2) - median(origin.unsaturated$D1) 	# Difference of sensitivity between D2 and D1.	
			}  
  		# Correction for the difference of sensitivity between D1 and D2
		if(origin > 0)  unsaturated$D2 <- 	unsaturated$D2 + origin	
  		if(origin < 0)  unsaturated$D1 <- 	unsaturated$D1 - origin	

  		aligned <- subset(unsaturated, D2 > (D1 - width * 10^4) & D2 < (D1 + width*10^4)) # filtering aligned particles (D1 = D2)
		focused <- subset(aligned, D2/fsc_small > (D1/fsc_small - notch) & D2/fsc_small < (D1/fsc_small + notch))# filtering focused particles (D1/fsc_small = D2/fsc_small)
  		filtered <- subset(focused, D1/fsc_small < notch | D2/fsc_small < notch) # filtering focused particles (D/fsc_small < notch)
  		filtered[,-c(1,2)] <-  10^((filtered[,-c(1,2)]/2^16)*3.5)
  		
  ########################
  
  if(do.plot){
	percent.opp <- round(100*nrow(filtered)/nrow(events),1)
	
	par(mfrow=c(2,2),pty='s')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
	
	events. <- events[1:round(nrow(filtered)),]
    
    plot(events.$D1, events.$D2, xlab="D1", ylab="D2", 
    		pch=20, cex=0.25, xlim=c(1,10^3.5), ylim=c(1,10^3.5), log='xy',
			col= densCols(log10(events.$D1), log10(events.$D2), colramp=.rainbow.cols),
			main="Filtering Aligned Particles")
	 mtext(paste("only", percent.opp,"% display"), side=3, cex=0.7,col='red')
	 abline(v=max(events[,"D1"]), h=max(events[,"D2"]), col='red',lwd=2)
	 abline(v=1, h=1, col='red',lwd=2)
		 par(new=T)
	 plot(1,1, xlim=c(0,2^16),ylim=c(0,2^16), bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, pch=NA)
	 abline(b=1, a=origin + width*10^4, col='red',lwd=2)
	 abline(b=1, a=origin -width*10^4, col='red',lwd=2)
	
	mtext(paste("Sensitivity Difference D2 =", round(origin,0)),outer=T,side=3, line=-1.5,font=2)
	mtext(paste("Width=", width),outer=T,side=3, line=-3,font=2)
	mtext(paste("Notch=", notch),outer=T,side=3, line=-4,font=2)

	aligned. <- subset(aligned, aligned$D1/aligned$fsc_small<2 & aligned$D2/aligned$fsc_small<2)[1:round(nrow(filtered)),]
    para1 <- aligned.$D1/aligned.$fsc_small
    para2 <- aligned.$D2/aligned.$fsc_small
    plot(x = para1,  y = para2, cex=.25, xlim=c(0,2), ylim=c(0,2),
         xlab='D1/fsc_small', ylab='D2/fsc_small', col = densCols(aligned.$D1/aligned.$fsc_small, aligned.$D2/aligned.$fsc_small,colramp=.rainbow.cols), pch=20, 
         main = 'Filtering Focused Particles',pty='s')    
	 mtext(paste("only", percent.opp,"% display"), side=3, cex=0.7,col='red')
     abline(v=notch, h=notch, col='red',lwd=2)
	 abline(b=1, a=notch, col='red', lwd=2)
	 abline(b=1, a=-notch, col='red', lwd=2)
		
   plot(filtered$fsc_small, filtered$pe, xlab="fsc_small", ylab="pe", 
    		pch=20, cex=0.25, xlim=c(1,10^3.5), ylim=c(1,10^3.5), log='xy',
			col= densCols(log10(filtered$fsc_small), log10(filtered$pe), colramp=.rainbow.cols),
			main=paste("Optimally Positioned Particles\n (",percent.opp,"% total events)"))
   plot(filtered$fsc_small, filtered$chl_small, xlab="fsc_small", ylab="chl_small", 
    		pch=20, cex=0.25, xlim=c(1,10^3.5), ylim=c(1,10^3.5), log='xy',
			col= densCols(log10(filtered$fsc_small), log10(filtered$chl_small), colramp=.rainbow.cols),
			main=paste("Optimally Positioned Particles\n (",percent.opp,"% total events)"))

  }
  
  
  return(filtered)
  
}



filterFile <- function(evt.path, width=1, notch=1, origin=NA, output.path=getCruisePath(evt.path)){

  path.pieces <- strsplit(evt.path,'/')
  year_day <- path.pieces[[1]][length(path.pieces[[1]])-1]
  julian_day <- substr(year_day,6,8)
  third_minute <- path.pieces[[1]][length(path.pieces[[1]])]

  message('loading')
  my.flow.frame <- readSeaflow(evt.path, transform=TRUE)


  out.evt.file.path <- .createOutputPath(evt.path, output.path)

  dir.create(dirname(out.evt.file.path))
  
  bitmap(paste(out.evt.file.path, '.opp.gif', sep=''), width=.plot.width, height=.plot.width)

  #def.par <- par(no.readonly = TRUE)
  #layout(rbind(c(1,2),c(1,3),c(4,5), c(4,5)))
  #par(mar=c(5,5,2,1),oma=c(1,1,1,1),pty='m')
  # plot(1,1,pch='',xlab='',ylab='',xaxt='n',yaxt='n',bty='n')  # TITLE of QUALITY CONTROL PLOt
  
  filter.frame <- filter(my.flow.frame, width=width, notch=notch, origin=origin, do.plot=TRUE)
  n <- nrow(filter.frame)
  message(paste(n, 'events found'))  

  message('writing')
  writeSeaflow(paste(out.evt.file.path,'.opp',sep=''), filter.frame)

    # ## now add a map plot so we can know where we are
    # latlong <- .getSDSlatlon(.getSDS(evt.path))
    # if(!is.null(latlong) & !any(sapply(latlong, is.na))){
      # cruise.track <- try(.makeSDSlatlonDF(read.delim(paste(output.path,'/sds.tab',sep=''))), silent=T) 
      # plotLatLongMap(latlong[1], latlong[2], cruise.track, margin=as.numeric(map.margin), pch=3, cex=2, lwd=2, col='red', main=paste("file: ",.getYearDay(evt.path),'/',basename(evt.path), sep=""))
      # }
   dev.off()	  


  #temporary tab output hack
  event_number <- 1:n

  if(FALSE){
    filter.frame2 <- cbind.data.frame(julian_day, third_minute, event_number, filter.frame[,-1])
    write.table(filter.frame2, paste(out.evt.file.path,'.opp.tab',sep=''),  quote=FALSE, sep='\t',row.names=FALSE, col.names=FALSE)
  }
}
