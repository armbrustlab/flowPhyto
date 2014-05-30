readSeaflow <- function(file.path, column.names = EVT.HEADER, column.size = 2, count.only=FALSE, transform=TRUE, add.yearday.file=FALSE){ 

  if(!(substr(file.path, nchar(file.path)-2, nchar(file.path)) %in% c('opp','evt')))
    warning("attempting to read a seaflow file that doesn't have an evt or opp extension")
  # reads a binary seaflow event file into memory as a dataframe
  if(!file.exists(file.path)){
    
    stop(paste("The file doesn't exist;", file.path))
    
  }else{
    
    ## initialize dimentional parameters
    n.bytes.header <- 4
    n.bytes.EOL <- 4
    n.columns <- length(column.names)
    n.extra.columns <- n.bytes.EOL / column.size  # number of 16 bit integers (2:10&0) in the EOL character
    n.int.columns <- n.columns + n.extra.columns
    n.bytes.file <- file.info(file.path)$size
    n.rows <- ((n.bytes.file - n.bytes.header) / column.size) / n.int.columns  
    n.events <- n.int.columns * n.rows

    ## open binary file for reading 
    con <- file(description = file.path, open="rb") 
    header <- readBin(con, integer(), n = 1, size = n.bytes.header, endian = "little")
    header.EOL <- readBin(con, integer(), n = 1, size = n.bytes.EOL, endian = "little")

    ## check number of events
    if(n.rows != header)
      stop(paste("ERROR: the predicted number of rows (",n.rows,") doesn't equal the header specified number of events (", header,")",sep=''))

    if(count.only){
     return(header) #return just the event count in the header
    }else{
      
      ## read the actual events
      integer.vector <- readBin(con, integer(), n = n.events, size = column.size, signed = FALSE, endian = "little")    
      ## reformat the vector into a matrix -> dataframe
      integer.matrix <- matrix(integer.vector[1:n.events], nrow = n.rows, ncol = n.int.columns, byrow=TRUE)
      integer.dataframe <- data.frame(integer.matrix[,1:n.columns])
      ## name the columns
      names(integer.dataframe) <- c(column.names)
      close(con) 
      
      ## Transform data to LOG scale
      
      if(transform) integer.dataframe[,EVT.HEADER[-c(1,2)]] <- 10^((integer.dataframe[,EVT.HEADER[-c(1,2)]]/2^16)*3.5)  
           
           
           
      if(add.yearday.file){
              integer.dataframe$file  <- getFileNumber(file.path)
              integer.dataframe$year_day <- .getYearDay(file.path)
      }
      
	
      return (integer.dataframe)
    }
  }
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


.writeFCSheader <- function(con, offsets){
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


writeSeaflow <- function(file.path, df, column.names = EVT.HEADER, unstransform=TRUE, type = "FCS"){
  if(!all(EVT.HEADER %in% names(df)))
    warning("attempting to read a seaflow file that doesn't have an evt or opp extension")

 
	## UNTRANSFORM LOG-SCALED DATA (BACK TO ORIGINAL DATA)
	if(unstransform) df[,EVT.HEADER[-c(1,2)]] <- (log10(df[,EVT.HEADER[-c(1,2)]])/3.5)*2^16

   ## open connection ##
  con <- file(description = file.path, open="wb")
 
  if(type == "FCS"){
    require("flowCore")

    ## CREATE A FLOWFRAME 
    x <- flowFrame(as.matrix(df))

 

    ## Grab information from SDS file
    sds <- read.delim(paste(dirname(dirname(file.path)),"/sds.tab",sep=""))
    file.number <- basename(sub(".evt.opp","",file.path))
    day <- basename(dirname(file.path))
    sds.line <- sds[sds$file == file.number & sds$day == day ,]

    time <- sds.line$time
    lat <- sds.line$lat
    long <- sds.line$long
    temp <- sds.line$temperature
    sal <- sds.line$salinity
    fluo <- sds.line$fluorescence
    red <- sds.line$bulk.red
    PAR <- sds.line$PAR
    flow.rate <- sds.line$flow.rate
    event.rate <- sds.line$event.rate

    # Total number of EVT detected
    total.event <- readSeaflow(sub('.opp','',file.path), count.only=TRUE)

    
    filename <- paste(day,"_",file.number,".fcs",sep="")

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
      
      description(x) <- list( '$CYT'="SeaFlow",
                '$INST'="University of Washington",
                '$LAB'="Armbrust Lab",
                '$FILE'=filename,
                '$TIME.UTC'= time,
                '$LATITUDE'=lat,
                '$LONGITUDE'=long,
                '$TEMP'=temp,
                '$SALINITY'=sal,
                '$CHLOROPHYLL'=fluo,
                '$BULK_RED'=red,
                '$PAR'=PAR,
                '$FLOW.RATE'=flow.rate,
                '$EVENT.RATE'=event.rate,
                '$ACQUISITION.TIME'=3,
                '$TOTAL.EVT'=total.event,
                '$FCSversion'="3.0")
                
      identifier(x) <- filename
      
      ld <- length(exprs(x)) * 2
      ctxt <- .collapseDesc(x)
      endTxt <- nchar(ctxt) + begTxt - 1
      endDat <- ld + endTxt
      endTxt <- endTxt + (nchar(endTxt + 1) - 1) + (nchar(endDat) - 1)
      endDat <- ld + endTxt
      description(x) <- list(`$BEGINDATA` = endTxt + 1, `$ENDDATA` = endTxt +  ld)
      ctxt <- .collapseDesc(x)
      offsets <- c(begTxt, endTxt, endTxt + 1, endTxt + ld, 0, 0)
      con <- file(paste(dirname(file.path), "/",filename,sep=""),open = "wb")
      on.exit(close(con))
      .writeFCSheader(con, offsets)
      writeChar(ctxt, con, eos = NULL)
      writeBin(as(t(exprs(x)), "integer"), con, size = 2, endian = 'little')
      writeChar("00000000", con, eos = NULL)
    print(paste("writing",filename))
    

  }


  if(type =="OPP"){
  ## write newline ##
   ## NEED TO ADD CHECK and REORDERING OF COLUMN NAMES
  n.bytes.header <- 4
  column.size <- 2
  EOL.double <- 10

  writeBin(as.integer(c(nrow(df),EOL.double)), con, size = n.bytes.header, endian = "little")

  ## write out end of line character
  #writeBin(10, con, size = 4, endian = "little") #done 2 lines above?

  ## construct a vector of integers from the dataframe with the EOL integers at the end of each line
  out.vect <- as.integer(unlist(t(cbind(df, EOL.double, 0))))
  
  out.vect <- out.vect[1:(length(out.vect)-2)] # hack to remove the last two \r\n characters (see below)

  ## write it out
  writeBin(out.vect, con, size = column.size)
  }
  close(con)
}





write.delim <- function(df, file, quote=FALSE, row.names=FALSE, sep='\t', ...){
	write.table(df, file,  quote=quote, row.names=row.names, sep=sep, ...)
}

