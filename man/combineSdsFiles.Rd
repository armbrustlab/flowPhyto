\name{combineSdsFiles}
\alias{combineSdsFiles}
\title{Combine the SDS files.}
\description{
The SDS files in each directory are concatenated into one dataframe via this function.
}
\usage{
combineSdsFiles(cruise.dir='.')
}
\arguments{
\item{cruise.dir}{Path to the cruise directory.}
}
\value{
A data frame representing the concatenation of all cruise subdirectory SDS files.
}
\examples{

seaflow.path <- system.file("extdata", 'seaflow_cruise', package="flowPhyto")
sds <- combineSdsFiles(seaflow.path)
sds
}

