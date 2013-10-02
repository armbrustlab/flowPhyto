library(flowPhyto, quietly=T)

cmdArgsToVariables() 

if(exists('evt.path')){ 

  filterFile(evt.path=evt.path, output.path=output.path, width=width, notch=notch, slope=slope, edge=edge)
  
}