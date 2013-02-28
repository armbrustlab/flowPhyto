\name{validatePopDef}
\alias{validatePopDef}
\title{Validate a Population Definition Dataframe.}
\description{
Validate the columns and values of the population definition dataframe passed to this function.
}
\usage{
validatePopDef(pop.def)
}
\arguments{
\item{pop.def}{Path to the raw event file to be filtered.}
}
\value{
a boolean indicating weither or not the pop def passed the validation check.
}
\examples{


seaflow.path <- system.file("extdata", 'seaflow_cruise', package="flowPhyto")
pop.def <- readPopDef(seaflow.path)
validatePopDef(pop.def)

}

