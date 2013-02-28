\name{census}
\alias{census}
\title{Cross tabulate a consensus vector}
\description{
Cross tabulate the population composition
}
\usage{
census(v, pop.def)
}
\arguments{
\item{v}{a consensus vector of population classifications.}
\item{pop.def}{A population (rows) definition dataframe with parameters (columns) for gating and clustersing.}
}
\value{
a one row, cross-tabulated dataframe of counts with one column for each population specified by the rows in the pop.def dataframe. Zeros are filled in for absent populations.
}
\examples{
opp.file.path <- system.file("extdata","seaflow_cruise","2011_001", "1.evt.opp", 
				package="flowPhyto")
pop.file.path <- system.file("extdata","seaflow_cruise","pop.def.tab", 
				package="flowPhyto")

opp <- readSeaflow(opp.file.path)
def <- readPopDef(pop.file.path)
pop <- classify(x=opp, pop.def= def)

census(v=pop$pop, pop.def=def)
}


