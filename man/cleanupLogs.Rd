\name{cleanupLogs}
\alias{cleanupLogs}
\title{Remove R Batch Output Files.}
\description{
Each step of the flowPhyto pipeline generates a large number of .Rout files that should be cleaned up.  This function does just that and optionally leaves the outputs with errors intact for troubleshooting purposes.
}
\usage{
cleanupLogs(log.dir='.', keep.erred=TRUE)
}
\arguments{
\item{log.dir}{The directory where the log files were written to.}
\item{keep.erred}{Whether or not to keep the log files with errors.}
}
\value{
none
}
\examples{


cleanupLogs('.')

}

