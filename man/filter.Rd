\name{filter}
\alias{filter}
\title{Filter Optimally Positioned Particles from an EVT dataframe}
\description{
The function normalizes the signals of the two position-sensitive detectors (D1 and D2) by the forward angle light scatter (fsc_small) signal to identify optimally positioned particles (OPP) from an EVT dataframe. Optionally the function outputs a control quality plot for OPP filtration.
}
\usage{
filter(events,  width=1, notch=1, slope=NULL, edge=1, do.plot=FALSE)
}
\arguments{
\item{events}{event dataframe}
\item{slope}{correction factor for the stream alignment. When stream is not properly aligned, aligned particles do not scatter light equally on D1 and D2 and therefore do not lie onto the 1:1 line on the scatter plot of D2 vs D1. By default, the value of the slope is calculated as the ratio of D2/D1}
\item{width}{the width of the gate to the sides of the 1:1 equal detector response defines the allowed error in particle trajectories across the width of the stream.}
\item{notch}{the correction factor for the sensitivity of FSC with respect to D1 and D2. Scattered light from focused particles is maximal at the forward scatter detector (FSC) and minimal at both position detectors. When the sensitivity of FSC and D1/D2 detectors is adjusted to respond equally to focused calibration particles, the FSC normalized by the signal of both position detectors must be lower than 1.} 
\item{edge}{location of the boundary layer between water/air. Particles located at the boundary layer scatter light that can be detected by the position detectors.}
\item{do.plot}{create a plot that showed the different steps for filtering out non-optimally positioned particles}

}
\value{
a optimal-position filtered event dataframe
}

\examples{
evt.file.path <- system.file("extdata","seaflow_cruise","2011_001", "1.evt", 
				package="flowPhyto")
evt <- readSeaflow(evt.file.path)
opp <- filter(evt)
summary(opp)
}


