\name{OPPtoFCSconversion}
\alias{OPPtoFCSconversion}
\title{Convert OPP files into FCS 3.0.}
\description{
This is a function that convert binary Optimally Positioned Particle (OPP) file(s) into a flowFrame and then save the flowFrame as a Flow Cytometry Standard (FCS) file(s) (3.0 version).
}
\usage{
OPPtoFCSconversion(opp.filelist, save.path = getCruisePath(opp.filelist),...)
}
\arguments{
\item{opp.filelist}{a list of one or multiple OPP files.}
\item{save.path}{path where fcs files will be saved.}
\item{...}{additional arguments passed to the flowFrame function (see flowCore package for details}

\value{
.FCS files
}

\examples{

opp.filelist <- system.file("extdata","seaflow_cruise","2011_001", "1.evt.opp", 
				package="flowPhyto")
				
y <- OPPtoFCSconversion(opp.filelist, save.path=NULL)
y

}

