\name{consensusFile}
\alias{consensusFile}
\title{Create a consensus from several classification vector files and cross tabulate the population composition}
\description{
Create a consensus from several classification vectors and cross tabulate the population composition
}
\usage{
consensusFile(opp.path, pattern='.[0-9]+-class.vct$', 
		output.path= paste(.createOutputPath(opp.path, 
			getCruisePath(opp.path)), ".consensus.vct", sep =""))
}
\arguments{
\item{opp.path}{Path to OPP event file.}
\item{pattern}{The suffix regular expression pattern used to find the n pass vector files for this opp file}
\item{output.path}{Path to the directory where you wish to output data.}
}
\value{
a one dimentional consensus vector file on disk and an invisible vector in memory
}
\examples{

cruise.nm <- 'seaflow_cruise'
temp.out.dir <- '.'

seaflow.path <- system.file("extdata", cruise.nm, package="flowPhyto")
file.copy(from=seaflow.path, to=temp.out.dir, recursive=TRUE)

opp.file.path <- system.file("extdata","seaflow_cruise","2011_001", "1.evt.opp", 
				package="flowPhyto")

consensusFile(opp.path=paste(temp.out.dir,'/',cruise.nm,'/2011_001/1.evt.opp',sep=''))
unlink(cruise.nm, recursive=TRUE)
}


