\name{plotCruiseStats}
\alias{plotCruiseStats}
\title{Plot a seaflow cruise statistics}
\description{
Read in the pop.stats.tab file and plot maps, or line plots of it and optionall sds info
}
\usage{
plotCruiseStats(cruise, x.var='map', y.vars=c('conc','chl_small','fsc_small'), 
		pops= c("ultra", "synecho"), sds.var=NULL,
		date.range=as.POSIXct(c("2009-01-01", "2099-12-31"), tz='UTC'), 
		output.path=paste(REPO.PATH, cruise,'/plots/',sep=''), ...)
}
\arguments{
\item{cruise}{Simplified cruise name (same name as the subdirectory in the seaflow data dir).}
\item{x.var}{X variable: Either map or lat, long or time.}
\item{y.vars}{Y variables: either conc, fluor, or fsc.}
\item{pops}{Which populattions to plot.  See the pop datastructure for abreviations to use.}
\item{sds.var}{Which of the sds variables to plot as a secondary axis in a line plot.}
\item{date.range}{date range.}
\item{output.path}{Path to the directory where you wish to output data.}
\item{...}{Additional parameters passed to plot.}
}
\value{
an overview statistics plot file is output to disk
}

\examples{

cruise <- system.file("extdata","seaflow_cruise", package="flowPhyto")

plotCruiseStats(cruise=cruise, output.path='.')

}


