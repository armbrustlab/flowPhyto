\name{summarizeFile}
\alias{summarizeFile}
\title{Compute aggregate statistics on a collection of opp or fcs files.}
\description{
perform aggregate statistics on a particular combination of filtered opp or fcs files for a particular population.
}
\usage{
summarizeFile(opp.paths, pop.names, output.path=getCruisePath(opp.paths[1]))
}
\arguments{
\item{opp.paths}{Path to the raw event file to be filtered.}
\item{pop.names}{Abreviated name of the population subset to be summarized.}
\item{output.path}{Path to the directory where you wish to output data.}
}
\value{
a single row dataframe (with header) file of per population aggregate statistics on both channels and log meta information
}
\examples{

example.cruise.name <- 'seaflow_cruise'
temp.out.dir <- '.'

seaflow.path <- system.file("extdata", example.cruise.name, package="flowPhyto")
file.copy(from=seaflow.path, to=temp.out.dir, recursive=TRUE)

opp.paths <- sapply(c(1,2,3), function(i) 
	system.file("extdata","seaflow_cruise","2011_001", paste(i,'.evt.opp',sep=''), 
			package="flowPhyto"))


summarizeFile(opp.paths, pop.names='nano', output.path='.')

## optionally create a resample list 
## to concatenate files conditionally on population concentrations
#resamp.list <- createResampleList(cruise.path, 
#			resample.min=500, resamp.concat.max=5)
unlink(example.cruise.name, recursive=TRUE)
}

