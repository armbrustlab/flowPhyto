library(flowPhyto, quietly=T)

cmdArgsToVariables() 

if(exists('opp.paths')){ 

 opp.paths.vect <- strsplit(opp.paths,',')[[1]]
  
  summarizeFile(opp.paths=opp.paths.vect, pop.names=strsplit(pop.names,',')[[1]], output.path=output.path)
}