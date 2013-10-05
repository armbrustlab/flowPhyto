.DecadeCorrection <-  function(x){
	t <- 10^((x/2^16)*3.5)
	return(t)
	}
	
	
	
.collapseDesc <- function (x){
    d <- description(x)
    d <- lapply(d, function(y) {
        if (length(y) == 0) 
            return(" ")
        else {
            if (is.matrix(y)) 
                return(y)
            else return(sub("^$", " ", y))
        }
    })
    d <- d[order(names(d))]
    spillName <- intersect(c("SPILL", "spillover"), names(d))
    if (length(spillName) > 0) {
        mat <- d[[spillName]]
        rNum <- as.character(nrow(mat))
        clNames <- paste(colnames(mat), sep = ",")
        vec <- paste(c(t(mat)), sep = ",", collapse = ",")
        d[spillName] <- paste(c(rNum, clNames, vec), sep = ",", 
            collapse = ",")
    }
    paste("/", iconv(paste(names(d), "/", sapply(d, paste, 
        collapse = " "), "/", collapse = "", sep = ""), to = "latin1", 
        sub = " "), sep = "")
}


.writeFCSheader <- function(con, offsets)
{
    seek(con, 0)
    writeChar("FCS3.0    ", con, eos=NULL)
    len <- length(offsets)/2
    for (i in seq_len(len)) {
        indx <- 2*(i-1) +1;
        val1 <- offsets[indx]
        val2 <- offsets[indx+1];
        st1 <- 8 - nchar(val1)
        st2 <- 8 - nchar(val2)
        if( nchar(val1) > 8 || nchar(val2) > 8){
             val1 <- val2 <- 0
        }
        writeChar(paste(paste(rep(" ", 8 - nchar(val1)), collapse=""), val1,
                        collapse="", sep=""), con, eos=NULL)
        writeChar(paste(paste(rep(" ", 8 - nchar(val2)), collapse=""), val2,
                        collapse="", sep=""), con, eos=NULL)
    }
     invisible()
}

OPPtoFCSconversion <- function(opp.filelist, save.path = getCruisePath(opp.filelist),...){

	
	cruise <- basename(getCruisePath(opp.filelist))

	system(paste("mkdir ", save.path,cruise,"/",sep=""))
	days <- unique(basename(dirname(opp.filelist)))
	for(d in days) system(paste("mkdir ", save.path, cruise,"/",d,sep=""))

	for(opp.file in opp.filelist){
		
		
		## DESCRIPTION
		year.day <- flowPhyto:::.getYearDay(opp.file)
		filename <- basename(opp.file)
		filenumber <- getFileNumber(opp.file)
		acquisition <- file.info(opp.file)$mtime
		last.modif <- file.info(opp.file)$ctime
		fcs.filename <- paste(cruise,".",year.day,".",filenumber,sep="")
		
		## READ OPP file
		opp <- readSeaflow(opp.file) 
		
		## GET RIDE OF PULSE WIDTH, D1 and D2	
		cleaned.opp <- opp[,c(1,5:10)] 
		
		## TRANFORM THE DATA 
		#cleaned.opp[,-1] <- .DecadeCorrection(cleaned.opp[,-1])

		## CREATE A FLOWFRAME	
		x <- flowFrame(as.matrix(cleaned.opp))
		
		## ANNOTATE		
		begTxt <- 58
		mk <- list(`$BEGINANALYSIS` = "0",
					`$BEGINDATA` = "0",
					`$BEGINSTEXT` = 0, 
        			`$BYTEORD` = "1,2,3,4",
        			`$DATATYPE` = "I",
        			`$ENDANALYSIS` = "0",
       				`$ENDDATA` = "0", 
      				`$ENDSTEXT` = "0",
        			`$MODE` = "L", 
        			`$NEXTDATA` = "0", 
        			`$PAR` = ncol(x), 
        			`$TOT` = nrow(x))
    	
    	pnb <- as.list(rep(16, ncol(x)))
    	names(pnb) <- sprintf("$P%sB", 1:ncol(x))
    	mk <- c(mk, pnb)
    
      	pne <- as.list(rep("4,0", ncol(x)))
    	names(pne) <- sprintf("$P%sE", 1:ncol(x))
    	mk <- c(mk, pne)
    
    	pnr <- as.list(c(unlist(list(range(x[,1])))[2], rep(2^16, ncol(x)-1)))
    	names(pnr) <- sprintf("$P%sR", 1:ncol(x))
    	mk <- c(mk, pnr)
    
    	pnn <- colnames(x)
   	 	names(pnn) <- sprintf("$P%sN", 1:ncol(x))
    	mk <- c(mk, pnn)
    
    	description(x) <- mk
    	
    	description(x) <- list( '$FCSversion'="3.0",
    							'$CYT'="SeaFlow",
 								'$CRUISE'=cruise, 
								'$YEARDAY'=year.day,
								'$DATE'= acquisition,
								'$LAST_MODIFIED'=last.modif,
								'$INST'="University of Washington",
								'$LAB'="Armbrust Lab",
								'$FIL'=fcs.filename)
								
    	identifier(x) <- fcs.filename
    	
    	ld <- length(exprs(x)) * 2
    	ctxt <- .collapseDesc(x)
    	endTxt <- nchar(ctxt) + begTxt - 1
    	endDat <- ld + endTxt
    	endTxt <- endTxt + (nchar(endTxt + 1) - 1) + (nchar(endDat) - 1)
  	 	endDat <- ld + endTxt
    	description(x) <- list(`$BEGINDATA` = endTxt + 1, `$ENDDATA` = endTxt +  ld)
    	ctxt <- .collapseDesc(x)
    	offsets <- c(begTxt, endTxt, endTxt + 1, endTxt + ld, 0, 0)
    	con <- file(paste(save.path,cruise,"/",year.day, "/", filenumber,".fcs",sep=""), open = "wb")
   		on.exit(close(con))
    	.writeFCSheader(con, offsets)
    	writeChar(ctxt, con, eos = NULL)
    	writeBin(as(t(exprs(x)), "integer"), con, size = 2, endian = 'little')
    	writeChar("00000000", con, eos = NULL)
		print(fcs.filename)
	}

return(x)
	
}