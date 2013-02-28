.DBcon <- function(){
  source("~/database.conf.R")
  for(var in c('.db.driver','.db.user','.db.pass','.db.name','.db.host',
               '.db.stats.tab.nm','.db.cruise.fkey.nm','.db.cruise.tab.nm','.db.sds.tab.nm'))
  assign(var, get(var), envir=.GlobalEnv)
  
  require(RPostgreSQL, quietly=T)

  dbConnect(dbDriver(.db.driver),
            user=.db.user, 
            password=.db.pass, 
            dbname=.db.name, 
            host=.db.host)
  
}

.loadCheck <- function(con, tab.nm, cruise.fkey.nm, cruise.id, update){
  ct <- dbGetQuery(con, paste("SELECT count(*) FROM", tab.nm, "WHERE", cruise.fkey.nm, "=",cruise.id))[[1]]

  if(ct > 0){
      if(update==FALSE){
        stop(paste("Found",ct,"prexisting", tab.nm, "records for this cruise. Try setting update to TRUE"))
      }else{
        dum <- dbGetQuery(con, paste("DELETE FROM", tab.nm, "WHERE", cruise.fkey.nm, "=",cruise.id))[[1]]
        print(paste("deleted", ct, tab.nm, "records from the database"))
      }
  }else{
    message(paste('Fresh', tab.nm, 'records for the database!'))
  }
}

.unqCnstrntCheck <- function(df, clmns){
  unqtgthr <- apply(df[,clmns], 1, function(r) gsub(' ','', paste(r, collapse='.')))
  nuniq.recs <- length(unique(unqtgthr))
  if(nuniq.recs != nrow(df)){
    unqtgthr.xtab <- table(unqtgthr)
    stop(paste("Unique constraint",paste(clmns,collapse='.'),"failed. \
    record count", nuniq.recs,"doesn't match number of rows",nrow(df),"of table. \
    Likely causes are the following records:", 
    paste(names(unqtgthr.xtab)[unqtgthr.xtab>1],collapse=', ')))
  }
}

.loadStats <- function(cruise, update=FALSE, con=.DBcon()){
  ## assumes you've already constructed a database table. with proper column names and types
  ## if not, please run the SQL script to create this table

  ## find the cruise id fromt the cruise table
  cruises <- dbReadTable(con, .db.cruise.tab.nm)
  cruise.id <- cruises$id[cruises$name_dir==cruise]

  .loadCheck(con, .db.stats.tab.nm, .db.cruise.fkey.nm, cruise.id, update)

  ## read in the actual table to import
  stats.tab <- read.delim(paste(REPO.PATH, cruise ,'/stats.tab',sep=''), stringsAsFactors=FALSE)

  .unqCnstrntCheck(df=stats.tab, clmns=c('day','resamp','pop'))

  # add something here to check for 2e+6 integers in evt
  exact.millions <- stats.tab$evt %% 1e+6
   if(sum(exact.millions==0) > 0)
     warning('you have several evt cts that are multiples of 1 million. This causes a known bug in postgres loading. specifically, the lines below:')
     print (1:nrow(stats.tab))[exact.millions==0]

  ## add the cruise id
  stats.tab.names <- names(stats.tab)
  stats.tab <- cbind.data.frame(stats.tab, as.integer(cruise.id))
  names(stats.tab) <- c(stats.tab.names, .db.cruise.fkey.nm)


  #add special db-only outlier identifier column
  stats.tab$outlier <- FALSE
  outliers.path <- paste(REPO.PATH, cruise ,'/outliers',sep='')
  if(file.exists(outliers.path)){
    outliers <- read.delim(outliers.path, stringsAsFactors=FALSE, header=FALSE)
    stats.tab[gsub(' ','',apply(stats.tab[,c('day','file')],1,paste, collapse=',')) %in% outliers[[1]],'outlier']<- TRUE
    stats.tab[gsub(' ','',paste(stats.tab[,  'day'        ],'*' ,         sep=',')) %in% outliers[[1]],'outlier']<- TRUE
  }
  
  message("loading statistics table to database")
  dbWriteTable2(con, .db.stats.tab.nm, stats.tab, pg.update.seq=TRUE, append=TRUE)
  
}	



.loadSDS <- function(cruise, update=FALSE,  con=.DBcon()){
  ## assumes you've already constructed a database table. with proper column names and types
  ## if not, please run the SQL script to create this table

  ## find the cruise id fromt the cruise table
  cruises <- dbReadTable(con, .db.cruise.tab.nm)
  cruise.id <- cruises$id[cruises$name_dir==cruise]

  ct <- dbGetQuery(con, paste("SELECT count(*) FROM", .db.sds.tab.nm, "WHERE", .db.cruise.fkey.nm, "=",cruise.id))[[1]]

  .loadCheck(con, .db.sds.tab.nm, .db.cruise.fkey.nm, cruise.id, update)

  ## read in the actual table to import 
  SDSs <- read.delim(paste(REPO.PATH, cruise ,'/sds.tab',sep=''), stringsAsFactors=FALSE)
  
  .unqCnstrntCheck(df=SDSs, clmns=c('day','file'))
  
  ## add the cruise id 
  sds.names <- names(SDSs)
  SDSs <- cbind.data.frame(SDSs, rep(as.integer(cruise.id),nrow(SDSs)))
  names(SDSs) <- c(sds.names, .db.cruise.fkey.nm)

  ## remove the null columns
  SDSs <- SDSs[,!(1:ncol(SDSs) %in% grep( 'NULL', names(SDSs)))] 

  names(SDSs) <- tolower(names(SDSs))

  message("loading sds table to database")
  dbWriteTable2(con, .db.sds.tab.nm, SDSs, pg.update.seq=TRUE, append=TRUE)

}	




.quotize <- function(var.vect)
  sapply(var.vect, function(x) paste("'",x,"'",sep=''))

  

.queryStats <- function(cruise, x.var='map',
                       y.vars = c('conc','fsc','chl'),
                       ranges = list(utc=c('2009-01-01','2030-01-01'),lat=c(-90,90),long=c(-180,180)),
                       populationsp = c('pico','ultra','nano','synecho','crypto')){
  if(FALSE){
    cruise <- 'Thompson'
    x.var <- 'time'
    y.vars <- c('conc','fsc','chl')
    ranges <- list(utc=c('2009-11-07','2009-11-09'),lat=c(-90,90),long=c(-180,180)) # .00,24.
    populations = c('pico','ultra','nano','synecho','crypto')
  }

  what <- paste("pop,", paste(y.vars,collapse=','))
  if(x.var == 'map'){
    what <- paste(what, ",lat,long")
  }else{
    what <- paste(what, " ,'", x.var,"'", sep='')
  }

  where <- paste(.db.stats.tab.nm,'.',.db.cruise.fkey.nm, '=',.db.cruise.tab.nm,".id AND name ='", cruise,"'", sep='')
  where <- paste(where, 'AND pop IN (', paste(.quotize(populations), collapse=','),')')
  for(range.name in names(ranges)){
    range <- ranges[range.name][[1]]
    if(all(sapply(range, is.character)))   #range.name == 'utc'
      range <- .quotize(range)
    where <- paste(where, 'AND', range.name, 'BETWEEN (', paste(range, collapse=' AND '),')')
  }

  
  sql <- paste("SELECT", what ,"FROM", .db.stats.tab.nm,',',.db.cruise.tab.nm, "WHERE", where,';')
  #rs <- dbSendQuery(con, sql)
  #df <- fetch(rs)
  #return(df)
  return(sql)

}

