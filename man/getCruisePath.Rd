\name{getCruisePath}
\alias{getCruisePath}
\title{Get the Cruise Directory.}
\description{
Retrieve the cruise directory path from a path string pointing to a file from that cruise
}
\usage{
getCruisePath(this.path, slash=TRUE)
}
\arguments{
\item{this.path}{Path to a file.}
\item{slash}{Boolean to indicate which if a slash should be used}
}
\value{
the cruise path's direcgtory
}
\examples{

path <- system.file("extdata","seaflow_cruise","2011_001","1.evt.opp", 
			package="flowPhyto")
getCruisePath(path)

}

