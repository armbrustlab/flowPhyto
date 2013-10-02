#######################
## CLASSIFY FUNCTION ##
#######################

classify <- function(x, pop.def=POP.DEF, varnames = CHANNEL.CLMNS.SM, numc=0, noise=0, auto.correction= FALSE, plot.cluster = FALSE, plot.assignment = FALSE, try.context='local',...){
  
  	    
    x$pop <- 0
    undefd <- subset(x, x$chl_small > as.numeric(noise))
    
	if(is.na(numc)){
		print("Calculating the optimal number of populations...")
		nc <- flowMeans(undefd, varNames = varnames, MaxN = NA , NumC = NA)
		NumC <- changepointDetection(nc@Mins)$MinIndex + 1
		print(paste("Clustering", NumC, "populations..."))
		res <- flowMeans(undefd, varNames = varnames, MaxN = NA, NumC = NumC)
				}else{
    
    if(as.numeric(numc) == 0){
	 	NumC <- nrow(pop.def)
		print(paste("Clustering", nrow(pop.def), "populations defined in pop.def table..."))
		res <- flowMeans(undefd, varNames = varnames, MaxN = NumC*2, NumC = NumC, ...) # CLUSTERING FUNCTION 
    	}
	 
	 if(as.numeric(numc) >= 1){
	 	NumC <- as.numeric(numc)
      	print(paste("Clustering", NumC, "populations..."))
      	res <- flowMeans(undefd, varNames = varnames, MaxN = NumC*2, NumC = NumC, ...) # CLUSTERING FUNCTION    
      		}
   		}
    
    ## Assigning each cell to a cluster number
    for(p in 1:NumC){
      df <- data.frame(subset(res@Label, res@Label == p))
      x[row.names(df),'pop'] <- as.numeric(p)
    }
    
    
    if(plot.cluster == TRUE){	
      ## Check Output of FlowMeans
      par(mfrow=c(3,2), pty='s')
      plot(undefd, res, c("fsc_small", "chl_small"), pch=1,xlim=c(0,2^16), ylim=c(0,2^16))
      plot(undefd, res, c("fsc_small", "pe"), pch=1,xlim=c(0,2^16), ylim=c(0,2^16))
      plot(undefd, res, c("chl_small", "pe"), pch=1,xlim=c(0,2^16), ylim=c(0,2^16))
      plot(undefd, res, c("fsc_small", "chl_big"), pch=1,xlim=c(0,2^16), ylim=c(0,2^16))
      plot(undefd, res, c("fsc_small", "fsc_perp"), pch=1,xlim=c(0,2^16), ylim=c(0,2^16))
      }
    
    
    ## Matching cluster number with cell population defined in pop.def.tab
    if(plot.assignment == TRUE){
      par(mfrow=c(ceiling(sqrt(NumC)),ceiling(sqrt(NumC))), pty="s")}
    
    for (p in 1:NumC){
      dat <- subset(x, pop == p)
      
      for(i in 1:nrow(pop.def)){
        xvar <- pop.def[i,"xvar"]
        yvar <- pop.def[i,"yvar"]
		xm <- median(dat[,xvar])
		ym <- median(dat[,yvar])
        
        if(!is.na(pop.def[i,"lim"])){
          if(xm > pop.def[i,"xmin"] & xm < pop.def[i,"xmax"] & ym > xm + pop.def[i,"lim"]){
            
            if(ym > pop.def[i,"ymin"] & ym < pop.def[i,"ymax"]){
              x[row.names(dat),'pop'] <- pop.def[i,"abrev"]
              
              if(plot.assignment == TRUE){
                plot(dat[,xvar], dat[,yvar], xlim=c(0,2^16), ylim=c(0,2^16), xlab=xvar, ylab=yvar, main=paste("cluster", p, "identified as", pop.def[i,"abrev"]));
                points(xm, ym, pch=3, col='red', lwd=2, cex=2);
                polygon(c(pop.def[i,"xmin"], pop.def[i,"xmax"], pop.def[i,"xmax"],pop.def[i,"xmin"]), c(pop.def[i,"ymin"],pop.def[i,"ymin"], pop.def[i,"ymax"], pop.def[i,"ymax"]), border="red", lwd=3);
                abline(b=1, a=pop.def[i,"lim"], col="red", lwd=3)}
              
              break
              }
          }
          
        }else{
          if(xm > pop.def[i,"xmin"] & xm < pop.def[i,"xmax"]){
            
            if(ym > pop.def[i,"ymin"] & ym < pop.def[i,"ymax"]){
              x[row.names(dat),'pop'] <- pop.def[i,"abrev"]
              
              if(plot.assignment == TRUE){
                plot(dat[,xvar], dat[,yvar], xlim=c(0,2^16), ylim=c(0,2^16), xlab=xvar, ylab=yvar, main=paste("cluster", p, "identified as", pop.def[i,"abrev"]));
                points(xm, ym, pch=3, col='red', lwd=2, cex=2);
                polygon(c(pop.def[i,"xmin"], pop.def[i,"xmax"], pop.def[i,"xmax"],pop.def[i,"xmin"]), c(pop.def[i,"ymin"],pop.def[i,"ymin"], pop.def[i,"ymax"], pop.def[i,"ymax"]), border="red", lwd=3)}
              
              break
            }
          }
        }	
      }
    }
    
	if(auto.correction == TRUE){    
	  	if(!any(x$pop == "beads") & any(pop.def$abrev == "beads")){
    			print("correct for Beads in any population")
    			xvar <- pop.def["beads", "xvar"]
    			yvar <- pop.def["beads", "yvar"]
    			beads1 <- subset(x, x[,yvar] > x[,xvar] + pop.def["beads", "lim"]) 
    			x[row.names(beads1),'pop'] <- "beads"
    			}	
    			
  	  	if(!any(x$pop == "synecho") & any(pop.def$abrev == "synecho")){
  	  		if(any(x$pop == "pico")){
  	  			print("correct for Synechococcus in Picoplankton population")
 	 	  		pico1 <- subset(x, x$pop == "pico")
  		  		xvar <- pop.def["synecho", "xvar"]
  		  		syn1 <- subset(pico1, pico1[,xvar] > pop.def["synecho", "xmin"])      
  		  		x[row.names(syn1),'pop'] <- "synecho"
  		  		}
  		  	}
  		  	
  		if(!any(x$pop == "crypto") & any(pop.def$abrev == "crypto")){
  			if(any(x$pop == "nano")){
  				print("correct for Cryptophytes in Nanoplankton population")
       		 	nano1 <- subset(x, x$pop == "nano")
       		 	xvar <- pop.def["crypto", "xvar"]
    			crypto1 <- subset(nano1, nano1[,xvar] > pop.def["crypto", "xmin"])
   		 		x[row.names(crypto1),'pop'] <- "crypto"
    			}
    		    		    		
    		if(any(x$pop == "beads") & any(pop.def$abrev == "crypto")){	   		
    			print("correct for Cryptophytes in Beads population")
    			beads1 <- subset(x, x$pop == "beads")
    			yvar <- pop.def["crypto", "yvar"]
    			xvar <- pop.def["crypto", "xvar"]
    			crypto2 <- subset(beads1, beads1[,yvar] > beads1[,xvar] + pop.def["crypto", "lim"])
    			x[row.names(crypto2),'pop'] <- "crypto"
    	 	  	}
    		}
    		  	
    	if(!any(x$pop == "nano") & any(pop.def$abrev == "nano")){
    		  	if(any(x$pop == "diatoms")){
    				print("correct for Nanoplankton in Pennate-like Diatom population")
    				diatoms1 <- subset(x, x$pop == "diatoms")
    				yvar <- pop.def["diatoms", "yvar"]
					xvar <- pop.def["diatoms", "xvar"]
    				nano1 <- subset(diatoms1, diatoms1[,yvar] < diatoms1[,xvar] + pop.def["diatoms", "lim"]) 
    				x[row.names(nano1),'pop'] <- "nano"
    				}
		
				if(any(x$pop == "lgdiatoms") & any(pop.def$abrev == "nano")){
					print("correct for Nanoplankton in large Pennate-like Diatom population")
					lgdiatoms1 <- subset(x, x$pop == "lgdiatoms") 
					yvar <- pop.def["lgdiatoms", "yvar"]
					xvar <- pop.def["lgdiatoms", "xvar"]
    				nano2 <- subset(lgdiatoms1, lgdiatoms1[,yvar] <  lgdiatoms1[,xvar] + pop.def["lgdiatoms", "lim"]) 
    				x[row.names(nano2),'pop'] <- "nano"
    				}
    	
    			if(any(x$pop == "ultra") & any(pop.def$abrev == "nano")){
					print("correct for Nanoplankton in Ultraplankton population")
					ultra1 <- subset(x, x$pop == "ultra")
					yvar <- pop.def["nano", "yvar"]
					xvar <- pop.def["nano", "xvar"]
					nano3 <- subset(ultra1, ultra1[,yvar] > -ultra1[,xvar] + 2*pop.def["nano","ymin"]) 
    				x[row.names(nano3),'pop'] <- "nano"
    				}
    			}	
  		  		
  		if(!any(x$pop == "ultra") & any(pop.def$abrev == "ultra")){
  			if(any(x$pop == "smdiatoms")){
			print("correct for Ultraplankton in small Pennate-like Diatom population")
			smdiatoms1 <- subset(x, x$pop == "smdiatoms") 
			yvar <- pop.def["smdiatoms", "yvar"]
			xvar <- pop.def["smdiatoms", "xvar"]
    		ultra1 <- subset(smdiatoms1, smdiatoms1[,yvar] < smdiatoms1[,xvar] + pop.def["smdiatoms", "lim"]) 
    		x[row.names(ultra1),'pop'] <- "ultra"
				}
  			}
  		
  		if(!any(x$pop == "cocco") & any(pop.def[,"abrev"] == "cocco")){
  				print("correct for Coccolithophores in any population")
				yvar <- pop.def["cocco", "yvar"]
				xvar <- pop.def["cocco", "xvar"]
				cocco1 <- subset(x, x[,yvar] > x[,xvar] + pop.def["cocco", "lim"] & x[,xvar] > pop.def["cocco", "xmin"] & x[,yvar] > pop.def["cocco", "ymin"])
				x[row.names(cocco1), 'pop'] <- "cocco"
				}	
		}		
			return(x)
		
}





####################################
## CLASSIFY FUNCTION for PIPELINE ##
####################################

classifyFile <- function(opp.path, concat.ct=3, output.path=getCruisePath(opp.path), def.path=paste(getCruisePath(opp.path),'pop.def.tab', sep=''),  varnames=CHANNEL.CLMNS.SM, numc=0, noise=0, auto.correction=FALSE){
  concat.ct <- as.integer(concat.ct)
  
  if(!is.na(def.path))
      pop.def <- readPopDef(def.path)
  else
      pop.def <- POP.DEF
      
  file.no <- getFileNumber(opp.path)
  yearday <- .getYearDay(opp.path)


  ## CONCATENATE FILES FOR CLUSTERING  
  opp.paths <- .getNextFiles(opp.path, n=concat.ct)  
  filtered <- do.call(rbind.data.frame, lapply(opp.paths, readSeaflow, add.yearday.file=TRUE))
   
  ## CLUSTERING 
  classified <- classify(x=filtered, pop.def=pop.def, try.context=paste(yearday,':',file.no), numc=numc, varnames=varnames, noise=noise, auto.correction = auto.correction)
  
 
  ## SPLIT CLASSIFICATION ($pop) VECTOR off evt.opp and back OUT INTO n NEW FILES
  for(path in opp.paths){
    class.path <- .createOutputPath(path, output.path)
    pop.vect <- subset(classified, file==getFileNumber(path), 'pop')
    write.delim(pop.vect, paste(class.path,'.',file.no,'-class.vct',sep=''))

  }
}



###################################
## Plot output Classify function ##
###################################

plotCytogram <- function(df, x.ax, y.ax, add.legend=FALSE, pop.def=POP.DEF, cex=0.5, pch=16, xlab=x.ax, ylab=y.ax, ...){
      
    plot(df[, x.ax], df[, y.ax], col='grey', pch=pch, cex=cex, xlim=c(0,2^16), ylim=c(0,2^16), xlab=xlab, ylab=ylab, ...)

    if(('pop' %in% names(df))){
      for(p in pop.def$abrev){
        df.p <- subset(df, pop==p)
        #print(paste('plotting population',p,':',nrow(df.p),'cells'))
        if(nrow(df.p)>2){
 #         try(points(df.p[,x.ax], df.p[, y.ax], col = densCols(df.p[,x.ax], df.p[,y.ax], col=colorRampPalette(c(hsv(h=rgb2hsv(col2rgb(pop.def[p,]$color))['h',], s=.05, v=.88), pop.def[p,]$color))), pch=pch, cex=cex))
        try(points(df.p[,x.ax], df.p[, y.ax], col = pop.def[p,]$color, pch=pch, cex=cex))
       }
      }
      if(add.legend){  
        legend('bottomright', legend=pop.def$title, pch=16, col=pop.def$color, bty='n')
      }
    }else{
      warning('there is no pop column to color populations (re: dataframe has not yet been classified)')
    }
}