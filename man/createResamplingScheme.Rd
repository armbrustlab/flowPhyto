\name{createResamplingScheme}
\alias{createResamplingScheme}
\title{Generate the Population Resampling scheme.}
\description{
.
}
\usage{
createResamplingScheme(cruise.path, resample.min=500, resamp.concat.max = 5)
}
\arguments{
\item{cruise.path}{Path cruise directory.}
\item{resample.min}{Minimum number of cells required for a population's resampling.}
\item{resamp.concat.max}{Maximum number of files allowed to be concatenated to get the resample minimum.}
}
\value{
a vector of comma delimted population named by the corresponding comma delimted year_day/file.
}
\examples{
seaflow.path <- system.file("extdata", 'seaflow_cruise', package="flowPhyto")
createResamplingScheme(seaflow.path)
}

