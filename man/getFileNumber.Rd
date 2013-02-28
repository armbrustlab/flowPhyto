\name{getFileNumber}
\alias{getFileNumber}
\title{Get the (original) integer file number of any seaflow repository file.}
\description{
Each seaflow EVT file is assigned a unique (per-directory) integer which get's carried on to subsequent processing steps.  This function extracts that number from any of the original or downstream files.
}
\usage{
getFileNumber(file.path)
}
\arguments{
\item{file.path}{Path to the file whose name you wish to extract a number from.}
}
\value{
an integer corresponding to the original event file number
}
\examples{

path <- system.file("extdata","seaflow_cruise","2011_001","1.evt.opp", 
			package="flowPhyto")
getFileNumber(path)

}

