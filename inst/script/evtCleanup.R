

library(flowPhyto, quietly=T)


library(caroline, quietly=T)

cmdArgsToVariables()

dir <- '/share/data/cruise/instrument/underway/seaflow/'

cruise <- 'CMOP_3'

files <- getCruiseFiles(paste(dir, cruise, sep=''))

for(file in files){

  ## read it
  evt <- readSeaflow(file)
  ## filter out false positive rows (all channels == 2^16/2

  bug.free <- apply(tmp, 1, function(x)  !all(x[c('D1','D2',CHANNEL.CLMNS)] == 2^15))
  
  if(table(bug.free)$TRUE !=0){
  out <- subset(tmp, bug.free)
  
  
  ## move old file
  
  ## write out new file
  
    ## UNFINISHED

  }
  
 }
