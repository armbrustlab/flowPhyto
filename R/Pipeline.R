
.cmd.rcmdbatch <- "R CMD BATCH --no-save --no-restore "
.complete.exts <- c('opp.gif','class.vct','.census.tab','stat.tab')  #used to check if jobs are done #NOT WORKING YET
.source.dir <- system.file('script',package='flowPhyto')

.removeFiles <- function(dir='.', prefix='', ext='\\.R\\.out', keep.erred=TRUE, keep.warned=FALSE){
  files <- list.files(dir, pattern=paste('^',prefix,'.*',ext,'$',sep='')) 
  for(file in files){
    path <- paste(dir,'/',file,sep='')
    if(keep.erred){
      vect <- readLines(path)
      if(length(grep('> proc\\.time()', vect[length(vect)-2]))<=0) # if not completed
        next
    }
    if(keep.warned){
      vect <- readLines(path)
      if(any(sapply(vect, function(x) length(grep('Warning message:', x))>0))) # if not warnings
        next
    }
    file.remove(path) 
  }
}

.countFinished <- function(paths, ext){
  files <- sapply(paste(paths,'*', ext, sep=''), Sys.glob)
  if(class(files)== 'list')
    files <- do.call(c,files)  
  return(length(files))
 }
      
.checkJobCompletion <- function(paths, ext, fold, type='', pct=.98){
    n.in.files <- length(paths)
    if(is.na(n.in.files))
      stop('no valid input files in directory.  check extensions')
    pct.cmplt <- 0

    while(pct.cmplt < pct){
      message(paste('Waiting for',type,'job steps to finish'))
      Sys.sleep(10) #* length(n.in.files))
      ct.fin <- .countFinished(paths, ext)
      ct.togo <- n.in.files * fold
      pct.cmplt <- ct.fin/ct.togo
      message(paste(ct.fin,'out of',ct.togo,'(',round(100*pct.cmplt,2),'%) of jobs have finished'))
    }
    return(ct.fin)
}


cleanupLogs <- function(log.dir='.', keep.erred=TRUE){
    ## clean up the cluster submission log & out files         
    .removeFiles(log.dir, prefix='',         ext='\\.R\\.out', keep.erred=keep.erred)
    .removeFiles('.',     prefix='STDIN\\.', ext='')
}

clearOutputs <- function(cruise.path='.', steps=1:4){
  # use this function to clean up the key progress indicator outputs of each step of the pipeline
  files <- unlist(sapply(.complete.exts[steps], function(e) getCruiseFiles(cruise.path,ext=e)))
  if(length(files)!= 0){
    tfvct<-sapply(files, file.remove)
    table(tfvct)
  }
}

.queueLimitCheck <- function(i, lim=5000, wait.time=300){
  if(i %% lim  == 0){
    message(paste('waiting',round(wait.time/60,1),'minutes for first', lim ,'jobs to finish'))
    Sys.sleep(wait.time)
    return(TRUE)
  }
  return(TRUE)
}

.processRawData <- function(cruise.dir='.', steps = 1:3, transform=TRUE, clust.concat.ct=3, filter.width=1.5, filter.notch=1, filter.edge=1, filter.slope=NA, classify.varnames=CHANNEL.CLMNS.SM, classify.func= 2, classify.numc=0, classify.noise=0, map.margin=2, output.path='.', log.dir='.', def.path=NULL, debug=FALSE, fresh=TRUE, cleanup=TRUE, pct=.97, parallel=FALSE, submit.cmd='qsub', range=NULL){
  
  
  if(length(steps)>3 | 4 %in% steps)
    stop('steps must be a vector of integers between 1 and 3')
  if(any(steps %in% c(1,2,3))){
    script.names <- c('Filter','Classify','Census')
    suffixes <- c('evt','opp','opp')
    complete.fold <- c(1, clust.concat.ct, 1)
    add.opts <- paste('--map.margin=', map.margin,' --concat.ct=', clust.concat.ct,' --def.path=', def.path,'--transform=', transform, ' --width=', filter.width,' --notch=', filter.notch, ' --slope=', filter.slope, ' --edge=', filter.edge, ' --func=',classify.func, ' --varnames=', paste(classify.varnames, collapse=','), ' --numc=', classify.numc, ' --noise=', classify.noise, sep='')
    

    
    ## farm per-event-file analysis out to a qsub compatible cluster scheduler
    for(i in steps){
    
      if(fresh)
        clearOutputs(cruise.dir, i)
    
      cmd.out <- paste('--output.path=',output.path, sep='')
      cmd.script <- paste(.source.dir,'/', script.names[i],'.R' ,sep='')
      cmd.log <- paste(log.dir,"/seaflow.", script.names[i], sep='')
      
      paths = getCruiseFiles(cruise.dir, ext=suffixes[i], range=range)
      for(j in seq(along=paths)){
        this.path <- paths[j]
  	cmd.filepath <- paste('--',suffixes[i],'.path=', this.path, sep='')
 	cmd.log.out <- paste(cmd.log, '.', .getYearDay(this.path),'.', getFileNumber(this.path), '.R.out', sep='')
	cmd <- paste(.cmd.rcmdbatch, cmd.filepath, cmd.out,  add.opts, cmd.script, cmd.log.out)

        if(parallel){
          cmd <- paste("echo '",cmd,"' | ", submit.cmd, sep='') 
          dummy <- .queueLimitCheck(j)        
        }
        
        message(paste('Job',j,'\n',cmd))
        Sys.sleep(0.1)
        system(cmd)

      }
      
      if(parallel)
        files.processed <- .checkJobCompletion(paths, .complete.exts[i], complete.fold[i], type=script.names[i], pct=pct)
      else
        files.processed <- 1
      if(script.names[i] == 'Census'){
        message('Concatenating Census tab files')
        write.delim(combineCensusFiles(cruise.dir),  paste(output.path,'/census.tab',sep=''))
      }
      if(cleanup){
        message(paste('Cleaning up',script.names[i],'log files'))
        cleanupLogs(log.dir)
      }    
    }
   
    return(files.processed)
  }else{
    return(1)
  }
}



.processStatistics <- function(cruise.path='.', resamp.vect, output.path='.', log.dir='.', cleanup=TRUE, pct=1.00, parallel=FALSE, submit.cmd='qsub'){ 
  cmd.out <- paste('--output.path=',output.path, sep='')


  outpaths <- vector()
  
  for(j in seq(along=resamp.vect)){
    resamp = resamp.vect[j]
    pops <- resamp
    paths <- paste(cruise.path,'/',strsplit(names(resamp),',')[[1]], '.evt.opp',sep='')
    
    cmd.pop <- paste(' --opp.paths=',paste(paths,collapse=','), ' --pop.names=',pops, sep='')
    cmd.log <- paste(log.dir,"/seaflow.",paste(sapply(paths,getFileNumber),collapse='+'),'.',pops,'.R.out', sep='')
    cmd.script <- paste(.source.dir,'/', 'Summarize.R', sep='')
    cmd <- paste(.cmd.rcmdbatch, cmd.pop, cmd.out, cmd.script , cmd.log)
    if(parallel){
      cmd <- paste("echo '",cmd,"' | ", submit.cmd, sep='')
      .queueLimitCheck(j)
    }
    message(paste('Job',j,'\n',cmd))
    system(cmd)
    outpaths <- c(outpaths, .createStatsPath(paths, pops))
    
  }

  if(parallel)
    files.processed <- .checkJobCompletion(outpaths,'', 1,type='Summarize', pct=pct)
  else
    files.processed <- 1
  message('Concatenating Statistics tab files')
  write.delim(.combinePopStatFiles(cruise.path),  paste(output.path,'/stats.tab',sep=''))
  
  if(cleanup){
    message('Cleaning up Statistics log files')
    cleanupLogs(log.dir)
  }
  
  return(files.processed)
}


pipeline  <- function(cruise.name='', repo=REPO.PATH, range=NULL, steps=1:4, pct=.97,  
	transform=TRUE, clust.concat.ct=3, resample.size=300, resamp.concat.max=10,
	filter.width=1.5, filter.notch=1, filter.edge=1, filter.slope=NA,
	classify.func=2, classify.varnames=CHANNEL.CLMNS.SM, classify.numc=0, classify.noise=0,
	map.margin=2, concat.sds=!is.na(match(1,steps)), load.to.db=FALSE,  preplot=FALSE, cleanup=TRUE,
	input.path=paste(repo, '/', cruise.name, sep=''),
	output.path=input.path,  log.dir=output.path,
	def.path=paste(input.path,'/', 'pop.def.tab',sep=''), parallel=TRUE, submit.cmd='qsub'){


  if(concat.sds){
    ## first collect the sds file and lump them into one file
    message("creating master sds")
    master.sds  <- combineSdsFiles(input.path)
    master.sds <- .cleanSDS(master.sds)
    write.delim(master.sds, paste(output.path,'/', 'sds.tab', sep=''))
  }

  
  ## filter/classify/census and resample events
  if(any(steps <= 4)){
    #setwd(output.path)
    prd <- .processRawData(input.path, steps[steps!=4], transform=transform, map.margin=map.margin,  filter.width=filter.width, filter.notch=filter.notch, filter.slope=filter.slope, filter.edge=filter.edge, classify.func=classify.func, classify.varnames=classify.varnames, classify.numc=classify.numc, classify.noise=classify.noise, clust.concat.ct=clust.concat.ct, output.path=output.path, log.dir=log.dir, def.path=def.path, debug=TRUE, cleanup=cleanup, pct=pct, parallel=parallel, submit.cmd=submit.cmd, range=range)
    
    
    if(4 %in% steps){
     resamp.vect <- createResamplingScheme(input.path, resample.size, resamp.concat.max)
     #sum(sapply(resamp.list, function(x) sum(!is.na(x))))                # 58407
     #sum(sapply(resamp.list, function(x) length(unique(x[!is.na(x)]))))  # 48054
     #length(unique( do.call(c,sapply(resamp.list, function(x) names(x[!is.na(x)])))))  #9565
     prd <- .processStatistics(input.path, resamp.vect, output.path=output.path, log.dir=log.dir, cleanup=cleanup, pct=pct, parallel=parallel, submit.cmd=submit.cmd)
    }
  }else{
    prd <- 1
  }
    
  if(prd <= 0){  # files.processed > 0
    stop(paste('Sorry, but only', prd,'files were processed. no finishing will take place'))
  }else{

    #message("creating master stats file")
    #master.stats <- .combinePopStatFiles(input.path, regexps=regexps)
    #write.delim(master.stats, paste(output.path,'/', 'pop.stats.tab', sep=''))
  
    if(load.to.db){
      message("Loading to database")

      con <- .DBcon()
      ## load the stats table to a/the database
      .loadStats(cruise.name, con=con, update=TRUE)
      ## load the SDS tables to a/the database
      .loadSDS(cruise.name, con=con, update=TRUE)

    }
  
    if(preplot){
      message("Preplotting Level 2")
      .prePlotLevel2(input.path, output.path, log.dir=paste(output.path,'/logs/',sep=''))
    }

  }
}
