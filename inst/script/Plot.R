
library(caroline, quietly=T)
sds.var = NULL
cmdArgsToVariables() 

if(exists('cruise')){ 
	library(flowPhyto, quietly=T)
	if(!all(sapply(flowPhyto:::.plot.args[1:3],exists)))
	    stop('not all of the plot ags variables exist')
	 
	  query.db <- FALSE

	  print(paste('--PARSING VARIABLES / FILTERING POP STATS TABLE (if applicable)--',date()))
	  y.vars.v <- parseArgString(y.vars, param.range=flowPhyto:::.outcome.vars, max.param.ct=length(flowPhyto:::.outcome.vars),min.param.ct=1)  
	  if(!(x.var %in% flowPhyto:::.indepnd.vars))
	    stop(paste("your x", x.var , "couldn't be found in the independent variable list (", paste(flowPhyto:::.indepnd.vars, collapse=','),")"))

	  pops.v <- parseArgString(populations, min.param.ct=1)  

	  date.range <- gsub('_',' ',date.range)
	  main.plot.title <- gsub('_',' ',title)


	  max.date.range <- c(as.POSIXct('2000-01-01', tz='UTC'), as.POSIXct('3000-01-01', tz='UTC'))
	  if(exists('date.range')){
	    date.range <- gsub('_',' ',date.range)
	    date.range.v <- parseArgString(date.range, param.range=max.date.range)
	  }else{
	    date.range.v <- max.date.range
	  }
  print(paste("plotCruiseStats(cruise='",cruise,"', x.var='",x.var,"', y.vars=c('",paste(y.vars.v,collapse='',''),"'), pops=c('",paste(pops.v,collapse='',''),"'), date.range=c('",paste(date.range.v,collapse='',''),"'), output.path='",output.path,"', sds.var='",sds.var,"')", sep=''))	  
 	  
  plotCruiseStats(cruise=cruise, x.var=x.var, y.vars=y.vars.v, pops=pops.v, date.range=date.range.v, output.path=output.path, sds.var=sds.var,  main=main.plot.title) #query.db=TRUE,
   
}
	