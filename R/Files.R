.pad <- function(vect,np){
 x <- as.numeric(vect)
 vect <- as.character(vect)
 for(p in 10^(1:np)){
   vect[x<p] <- paste("0",vect[x<p],sep="")
 }
 vect
}

.getYearDay <- function(file.path){
  
  dir.path.vect <- strsplit(dirname(file.path), '/')[[1]]
  if(length(dir.path.vect)==1) #for the sds file column (eg. 'sds_2010_128_209')
    return(paste(strsplit(file.path, '_')[[1]][c(2,3)], collapse='_'))
  else  # grabbing it from the full file path (eg './Thompson/2010_128/209.evt')
    return(dir.path.vect[length(dir.path.vect)])

}

.setYearDay <- function(file.path, year=NA, day=NA, add=FALSE){
  dir.path.vect <- strsplit(file.path, '/')[[1]]
  
  file.name <- dir.path.vect[length(dir.path.vect)]
  old.year.day <- dir.path.vect[length(dir.path.vect)-1]
  if(length(dir.path.vect) == 2)
    pre.path = ''
  else
    pre.path <- paste(paste(dir.path.vect[1:(length(dir.path.vect)-2)],collapse='/'),'/',sep='')
  
  old.y.d.v <- strsplit(old.year.day,'_')[[1]]
  old.year <- as.integer(old.y.d.v[1])
  old.day <- as.integer(old.y.d.v[2])
  

  if(is.na(day)){
     day = old.day
  }else{
     if(add){
     	if(old.day == 365){
     	  day = 1
     	  year = 1     	  
     	}else{
          day = old.day + day
        }
      }
  }
  
  if(is.na(year)){
       year = old.year
    }else{
       if(add)
          year = old.year + year
  }
  
  if(day > 365 | day < 1)
    stop('day is out of 1 to 365 julian day range')
  new.year.day <- paste(year, '_', .pad(day,2), sep='')
  return(paste(pre.path, new.year.day,'/', file.name, sep=''))
  
}

getFileNumber <- function(file.path){

  if(length(grep('sds_[0-9]{4}_[0-9]{1,3}_[0-9]{1,9}',file.path))==1){ #for the sds file column (eg. 'sds_2010_128_209')
    return(as.integer(strsplit(file.path, '_')[[1]][4]))
  }else{
    file.name <- basename(file.path)
    re.num <- regexpr('^[0-9]{1,9}', file.name)
    file.number <- substr(file.name, re.num, attr(re.num, 'match.length'))
    ##file.number <- sub('.evt.*','', file.name)
    return(as.integer(file.number))
  }
}


.setFileNumber <- function(file.path, file.number){
  
  dir <- dirname(file.path)
  file.name <- basename(file.path)
  
  number <- as.integer(file.number)
  new.file.name <- sub('[0-9]{1,9}', number, file.name)
  
  return(paste(dir,'/', new.file.name, sep=''))
  
}

.getDirMaxFileNumber <- function(file.path){
  ext <- tail(strsplit(file.path,'\\.')[[1]],1)
  pat <- paste(ext,'$',sep='')
  num <- max(getFileNumber(list.files(dirname(file.path), pattern=pat)))
  return(num)
}


.getNextFiles <- function(file.path, n=3, gap=2){
  files <- file.path
  begin <- getFileNumber(file.path)
  i <- 0 # file iterator
  jumps <- 0 # gap jumping iterator
  while(i < n+jumps  &  length(files) < n  &  jumps <= gap){
    print(i)
    i <- i + 1    
    next.file <- .setFileNumber(file.path, begin+i)
    
    if(file.exists(next.file)){
      files <- c(files, next.file)
      jumps <- 0 # reset the jumps
    }else{
      if(begin+i >= .getDirMaxFileNumber(next.file) ){
        ## end of this dir. look in next days' dir
        next.day.file.path <- .setFileNumber(.setYearDay(file.path, day=1, add=TRUE), 1)
        if(file.exists(next.day.file.path))
          files <- c(files, .getNextFiles(next.day.file.path,n-i))
        i <- n+gap # end the loop
      }else{
        ## missing file. jump it!
        jumps <- jumps + 1
      }
      
    }
  }
  numbers <- sapply(files, getFileNumber)
  names(files) <- numbers
  return(files)
}


getCruisePath <- function(this.path, slash=TRUE){
  dir.path.parts <- strsplit(dirname(this.path),"/")[[1]]
  out <- paste(dir.path.parts[-length(dir.path.parts)], collapse='/')
  if(out =='')
    out <- '.'
  if(slash)
    out <- paste(out,'/',sep='')

  return(out)
}


.createOutputPath <- function(in.path, out.dir)
    sub(getCruisePath(in.path), paste(out.dir,'/',sep=''), in.path, fixed = TRUE)
    


getCruiseFiles <- function(cruise.dir='.', prefix='[0-9]{1,9}', ext='\\.evt', range=NULL){
    paths = vector()
    byday.subdirs <- list.files(cruise.dir, pattern="^[0-9]{4}_[0-9]{3}$") 
    if(length(byday.subdirs)==0)
      stop(paste('no directories found in ', cruise.dir))

    if(!is.null(range)){
      if(is.null(names(range)) | length(range)!=2 | !is.numeric(range))
        stop('range parameter should be a (year_day) named (file number) integer vector of length = 2')
      byday.subdirs <- byday.subdirs[byday.subdirs >= names(range)[1] & byday.subdirs <= names(range)[2]]
    }
  
    for(day.dir in byday.subdirs){
      day.path <- paste(cruise.dir,'/' ,day.dir, sep='')
      files <- list.files(day.path, pattern=paste('^',prefix,'.*',ext,'$',sep=''))

      files <- files[order(getFileNumber(files))] # reorder properly
      finos <- getFileNumber(files)
                     
      if(!is.null(range) & day.dir %in% names(range)){
        if(day.dir == names(range)[1]){
          files <- files[finos >= range[1]]
          finos <- getFileNumber(files) # recalculate because vector size changes
        }  
        if(day.dir == names(range)[2])
          files <- files[finos <= range[2]]        
      }
      
      if(length(files)==0)
        warning(paste('no files found in directory', day.dir))
      for(this.file in files){
	paths <- c(paths, this.path <- paste(day.path,'/', this.file, sep=''))

      }
    }
    return(paths)
}


.getYear <- function(year_day)
  as.integer(sapply(year_day, function(s) strsplit(s,'_',)[[1]][1]))
    
.getDay <- function(year_day)
  as.integer(sapply(year_day, function(s) strsplit(s,'_',)[[1]][2]))


      #days <- .getDay(byday.subdirs)
      #range.d <- .getDay(names(range))
