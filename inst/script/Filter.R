
library(caroline, quietly=T)
cmdArgsToVariables() 

if(exists('evt.path')){ 

  library(flowPhyto, quietly=T)
  
  filterFile(evt.path=evt.path, output.path=output.path, width=width, notch=notch, slope=slope, edge=edge)
  
}