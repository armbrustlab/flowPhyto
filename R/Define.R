

## DEFINE DEFAULT POPULATION DEFINITION PARAMETERS ##
POP.DEF <- data.frame(
                  abrev = I(c('beads','synecho','crypto','cocco','diatoms','prochloro','ultra','nano','pico')),
                  title = I(c('Beads','Synechococcus','Cryptophyte','Coccolithophores','Elongated','Prochlorococcus','Ultraplankton','Nanoplankton','Picoplankton')),
                  xmin  = c(0.5, 0.5, 3.5, 2.0, 2.0, 1.0, 2.0, 3.5, 1.0) * 10^4,
                  ymin  = c(2.0, 0.5, 3.5, 2.0, 3.0, 1.0, 2.0, 3.0, 1.5) * 10^4,
                  xmax  = c(6.5, 3.5, 6.5, 4.5, 6.5, 2.0, 4.5, 6.5, 3.0) * 10^4,
                  ymax  = c(6.5, 3.0, 6.5, 4.5, 6.5, 2.0, 4.5, 6.5, 3.5) * 10^4,
                  color = I(c('gray40','tan2','tomato3','blueviolet','gold','violetred4','palegreen3','darkcyan','lightseagreen')),
                  xvar  = I(c('chl_small', rep('pe', 2),rep('fsc_small', 6))),
                  yvar  = I(c('pe',rep('chl_small', 2), 'fsc_perp', 'chl_big', rep('chl_small', 2), 'chl_big', 'chl_small')),
                  u.co  = c(0.05, 0.25, 0.25, 0.5, 0.25, 0.50, 0.50, 0.50, 0.50),
                  lim    = c(0.5, 0.5, 0.5, 0.5, 0.5, NA, NA, NA, NA) * 10^4
                  )

row.names(POP.DEF) <- as.character(POP.DEF$abrev)



readPopDef <- function(pop.def.tab.path){
      ## check to see if there is an externally defined pop definition table 
      if(file.info(pop.def.tab.path)$isdir)
          pop.def.tab.path <- paste(pop.def.tab.path,'/pop.def.tab',sep='') 
      if(file.exists(pop.def.tab.path)){
        pop.def <- read.delim(pop.def.tab.path, as.is=TRUE)
        if(!validatePopDef(pop.def))
            stop('This is not a valid pop.def file.  please read the documentation conderning proper format')
        rownames(pop.def) <- as.character(pop.def$abrev)
        return(pop.def)
      }else{ #if there is not one, write the hard coded one above
        warning('No pop.def.tab file found. Writing hardcoded one into specified directory')
        write.table(POP.DEF, pop.def.tab.path, quote=FALSE, sep='\t', row.names=TRUE)
        return(POP.DEF)
    }
}


validatePopDef <- function(pop.def){

    valid <- TRUE
    if(!all(names(POP.DEF) %in% names(pop.def))){
        warning("not all of the names in the default pop def match those in your custom one")
    	valid <- FALSE
    }
    if(!all(c(levels(pop.def$xvar),levels(pop.def$yvar)) %in% CHANNEL.CLMNS)){
        warning("not all of the levels of your x & y var colums of pop def match the global CHANNEL.CLMNS")
        valid <- FALSE
    }	
    
        
    return(valid)
}