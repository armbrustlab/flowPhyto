library(flowPhyto, quietly=T)

cmdArgsToVariables() 

if(exists('opp.path')){ 

  varnames <- parseArgString(varnames, max.param.ct=32, param.range=CHANNEL.CLMNS)
  
  classifyFile(opp.path=opp.path, concat.ct=concat.ct, output.path=output.path, def.path=def.path, varnames=varnames, numc=numc, noise=noise)

}
