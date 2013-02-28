\name{joinSDS}
\alias{joinSDS}
\title{join this dataframe to the corresponding SDS log file entry.}
\description{
perform aggregate statistics on a particular combination of filtered opp or fcs files for a particular population.
}
\usage{
joinSDS(x, opp.paths)
}
\arguments{
\item{x}{a OPP dataframe.}
\item{opp.paths}{the named paths of the SDS files (names should be character strings of the corresponding file)}
}
\value{
a dataframe of the merge between the events and the SDS log info.
}
\examples{


opp.path <- system.file("extdata", 'seaflow_cruise','2011_001','1.evt.opp', package="flowPhyto")
opp <- readSeaflow(opp.path, add.yearday.file=TRUE)

opp.join <- joinSDS(opp, caroline::nv(opp.path, 1))
summary(opp.join)
}

