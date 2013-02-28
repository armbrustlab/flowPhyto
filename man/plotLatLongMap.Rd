\name{plotLatLongMap}
\alias{plotLatLongMap}
\title{plot file location on a map}
\description{
plot the location of the file (longitude and latitude) on a map
}
\usage{
plotLatLongMap(lat, long, track=NULL, margin=0.1, col='red', legend=NULL, pch=20, cex=1.5, lwd=1, lty=2, xlim=NULL, ylim=NULL, zlab=NA, xlab="Longitude (deg W)",ylab="Latitude (deg N)",...)
}
\arguments{
\item{lat}{vector of latitudes.}
\item{long}{vecor of longitudes.}
\item{track}{a dataframe of the cruise track (from the sds file) with longitude (long) and latitude (lat) named accordingly.}
\item{margin}{margin of longitude and latitude around edges of current position.}
\item{col}{vector of colors of the cruise track.}
\item{legend}{vector of named colors for use in a heatmap color legend.}
\item{pch}{point character}
\item{cex}{character expansion}
\item{lwd}{line width}
\item{lty}{line type}
\item{xlim}{limit of x axis}
\item{ylim}{limit of y axis}
\item{zlab}{the label for the values of the heatmap color gradient passed as 'col'}
\item{xlab}{the x-axis label}
\item{ylab}{the y-axis label}
\item{...}{additional arguments to be passed to the plot function}
}

\examples{

stat.tab <- system.file("extdata","seaflow_cruise","stats.tab", 
				package="flowPhyto")
stats <- read.delim(stat.tab)
plotLatLongMap(stats$lat[10], stats$long[10], track=stats)

}


