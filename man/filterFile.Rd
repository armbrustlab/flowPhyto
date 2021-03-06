\name{filterFile}
\alias{filterFile}
\title{Filter an EVT file and output an Optimally Positioned Particles (OPP) file}
\description{
The function normalizes the signals of the two position-sensitive detectors (D1 and D2) by the forward angle light scatter (fsc_small) signal to identify optimally positioned particles (OPP) from an EVT file
}
\usage{
filterFile(evt.path, width=1, notch=1, origin=NA, output.path=getCruisePath(evt.path))

}
\arguments{
\item{evt.path}{path to the raw EVT file to be filtered.}
\item{origin}{correction factor for the stream alignment. When stream is not properly aligned, aligned particles do not scatter light equally on D1 and D2 and must be corrected. By default (NA), the value of the origin is calculated as the difference between the median value of D2 with respect to the median value of D1}
\item{width}{the width of the gate to the sides of the 1:1 equal detector response defines the allowed error in particle trajectories across the width of the stream.}
\item{notch}{the correction factor for the sensitivity of FSC with respect to D1 and D2. Scattered light from focused particles is maximal at the forward scatter detector (FSC) and minimal at both position detectors. When the sensitivity of FSC and D1/D2 detectors is adjusted to respond equally to focused calibration particles, the FSC normalized by the signal of both position detectors must be lower than 1.} 
\item{output.path}{path to the directory where you wish to output OPP file.}
}

\value{
a seaflow opp evt file and a plot of the filtration process
}
\examples{

example.cruise.name <- 'seaflow_cruise'
temp.dir <- '.'
seaflow.path <- system.file("extdata", example.cruise.name, package="flowPhyto")
file.copy(from=seaflow.path, to=temp.dir, recursive=TRUE)

filterFile(paste(temp.dir,'/',example.cruise.name,'/2011_001/1.evt',sep=''))
unlink(example.cruise.name, recursive=TRUE)
}


