\name{readPopDef}
\alias{readPopDef}
\title{Read the Population Definition File.}
\description{
Read the population definition file into memory from disk
}
\usage{
readPopDef(pop.def.tab.path)
}
\arguments{
\item{pop.def.tab.path}{Path to the population definition file or the cruise directory.}
}
\value{
a dataframe of population definition parameters
}
\examples{

seaflow.path <- system.file("extdata", 'seaflow_cruise', package="flowPhyto")
readPopDef(seaflow.path)

}
\seealso{
POP.DEF
}
