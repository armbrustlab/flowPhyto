\name{censusFile}
\alias{censusFile}
\title{Create a consensus from several classification vector files and cross tabulate the population composition}
\description{
Create a consensus from several classification vectors and cross tabulate the population composition
}
\usage{
censusFile(opp.path, map.margin=2, output.path=getCruisePath(opp.path), 
		def.path=paste(getCruisePath(opp.path),'pop.def.tab',sep=''))
}
\arguments{
\item{opp.path}{Path to OPP event file.}
\item{map.margin}{Margin in latitude/longitude around the map plots.}
\item{output.path}{Path to the directory where you wish to output data.}
\item{def.path}{Path to the file that defines how to gate & cluster the events into populations.}
}
\value{
a one dimentional consensus vector file (one column) 
and a one dimentional cross tabulation file (one row)
both writen to disk
}

\examples{
example.cruise.name <- 'seaflow_cruise'
temp.dir <- '.'

seaflow.path <- system.file("extdata", example.cruise.name, package="flowPhyto")
file.copy(from=seaflow.path, to=temp.dir, recursive=TRUE)

#opp.in.path <- system.file("extdata","seaflow_cruise","2011_001", "1.evt.opp", 
#					package="flowPhyto")

opp.out.path <- paste(temp.dir,'/',example.cruise.name,'/2011_001/1.evt.opp',sep='')
censusFile(opp.path=opp.out.path, map.margin=.5)
unlink(example.cruise.name, recursive=TRUE)
}


