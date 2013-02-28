\name{plotCytogram}
\alias{plotCytogram}
\title{Plot a Phytoplankton Cytogram}
\description{
Plot a Phytoplankton Cytogram
}
\usage{
plotCytogram(df, x.ax, y.ax, add.legend=FALSE, pop.def=POP.DEF, cex=0.5, pch=1, xlab = x.ax, ylab = y.ax, ...)
}

\arguments{
\item{df}{a dataframe of events (rows) and channels (columns).}
\item{x.ax}{column to plot in the x.axis}
\item{y.ax}{column to plot in the y.axis}
\item{add.legend}{should the plot automatically generate a legend}
\item{pop.def}{A population (rows) definition dataframe with parameters (columns) for gating and clustersing.}
\item{cex}{character expansion for the points.  Undefined background points are 1/3rd of the foreground points.}
\item{pch}{point character}
\item{xlab}{label for the x axis. by default equals x.ax}
\item{ylab}{label for the y axis. by default equals y.ax}
\item{...}{other paramters passed to plot}
}
\value{
a cytogram plot
}

\examples{
opp.file.path <- system.file("extdata","seaflow_cruise","2011_001", "1.evt.opp", 
			package="flowPhyto")
pop.file.path <- system.file("extdata","seaflow_cruise","pop.def.tab", 
			package="flowPhyto")

opp <- readSeaflow(opp.file.path)
def <- readPopDef(pop.file.path)
pop <- classify(x=opp, pop.def= def)

# Visualize the result of Classify using the function plotCytogram()
plotCytogram(pop, "fsc_small","chl_small", pop.def= def)
}


