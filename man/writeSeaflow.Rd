\name{writeSeaflow}
\alias{writeSeaflow}
\title{Write A SeaFlow File}
\description{
writes a binary seaflow event file from a dataframe in memory
}
\usage{
writeSeaflow(file.path, df, column.names = EVT.HEADER, linearize=TRUE)
}
\arguments{
\item{file.path}{System path to the binary seaflow event file.}
\item{df}{SeaFlow dataframe in memory to be written to disk.}
\item{column.names}{Names of the columns. By default it uses the global variable 'EVT.HEADER'}
\item{linearize}{Linearize log transformed data (use only if data were previously log transformed, see readSeaflow( ) function for more details).}
}

\examples{
opp.path <- system.file("extdata","seaflow_cruise","2011_001", "1.evt.opp",
				package="flowPhyto")
opp <- readSeaflow(opp.path)
writeSeaflow('./tmp.seaflow', opp)
Sys.sleep(30) # for windows build
}


