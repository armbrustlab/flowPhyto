\name{getCruiseFiles}
\alias{getCruiseFiles}
\title{Get all of the files from a cruise.}
\description{
This is a convenience function to grab all of the raw data files from the all of the julian day sub directories in a seaflow cruise directory.
}
\usage{
getCruiseFiles(cruise.dir='.', prefix='[0-9]{1,9}', ext='\\\\.evt', range=NULL)
}
\arguments{
\item{cruise.dir}{Path to cruise.}
\item{prefix}{a prefix to add to the files you wish to list.}
\item{ext}{extension of the files of interest.}
\item{range}{A named, two-integer vector specifying the start and end (inclusive) range for subsetting the input files used in each analysis step (with the exception of summarize). Values should be a (evt/opp) file numbers and names should be strings corresponding to the year_julianday directory names.  The nv() function is useful for creating this vector.}
}
\value{
a vector of cruise file names
}
\examples{

path <- system.file("extdata","seaflow_cruise", package="flowPhyto")
getCruiseFiles(path)

}

