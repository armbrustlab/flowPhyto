\name{readConsensusFile}
\alias{readConsensusFile}
\title{Read a Consensus File.}
\description{
Read a population classification consensus vector file from disk.
}
\usage{
readConsensusFile(path)
}
\arguments{
\item{path}{path to the consensus file}
}
\value{
a concensus vector in memory
}
\examples{

opp.path <-  system.file("extdata","seaflow_cruise","2011_001", '1.evt.opp',
			package="flowPhyto")

v <- readConsensusFile(opp.path)
table(v)
}

