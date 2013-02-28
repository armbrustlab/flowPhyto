\name{combineCensusFiles}
\alias{combineCensusFiles}
\title{Combine the Census Files.}
\description{
The census step generates a single row of population counts per EVT file.  This function collects and concatenates these into a single dataframe.
}
\usage{
combineCensusFiles(cruise.dir='.')
}
\arguments{
\item{cruise.dir}{Path the cruise directory.}
}
\value{
a dataframe of counts per EVT file (rows) and per population (columns) 
}
\examples{

seaflow.path <- system.file("extdata", 'seaflow_cruise', package="flowPhyto")
census <- combineCensusFiles(seaflow.path)
census

}

