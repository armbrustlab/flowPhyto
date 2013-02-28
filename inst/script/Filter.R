
library(caroline, quietly=T)
cmdArgsToVariables() 

if(exists('evt.path')){ 

  library(flowPhyto, quietly=T)
  
  filterFile(evt.path=evt.path, output.path=output.path, map.margin=map.margin, width=width,  notch=notch, edge=edge)
  
}