\name{summarize}
\alias{summarize}
\title{Summarize phytoplankon event & log parameters}
\description{
perform aggregate statistics on a particular combination of events (in a dataframe) per population.
}
\usage{
summarize(x, channel.clmns = CHANNEL.CLMNS, opp.paths.str='1,2,3')

}
\arguments{
\item{x}{event dataframe}
\item{channel.clmns}{flow cyometry event-level chanels columns in the dataframe (x) on which to perform aggreate statistics (mean and standard deviation).}
\item{opp.paths.str}{a comma delimited string of file names}
}
\value{
returns an aggregate statistics dataframe
}

\examples{

opp.paths <- sapply(c(1,2,3), function(i) 
      system.file("extdata","seaflow_cruise","2011_001", paste(i,'.evt.opp',sep=''), 
          package="flowPhyto"))

filter.df <- do.call(rbind.data.frame, lapply(opp.paths, readSeaflow, 
						add.yearday.file=TRUE))
consen.df <- do.call(rbind.data.frame, lapply(opp.paths, readConsensusFile)) 
classed <- cbind.data.frame(filter.df, consen.df)


nrow.opp <- sapply(opp.paths, function(p) readSeaflow(              p , count.only=TRUE))
nrow.evt <- sapply(opp.paths, function(p) readSeaflow(sub('.opp','',p), count.only=TRUE))
classed$opp <- rep(nrow.opp, times=nrow.opp)
classed$evt <- rep(nrow.evt, times=nrow.opp)
  
summary(classed)

class.jn <- do.call(rbind.data.frame, by(classed, list(classed$file), 
			joinSDS, caroline::nv(opp.paths, sapply(opp.paths, getFileNumber)) ))
summarize(class.jn,  opp.paths.str=paste(names(opp.paths), collapse=','))

}


