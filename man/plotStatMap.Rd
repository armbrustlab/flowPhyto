\name{plotStatMap}
\alias{plotStatMap}
\title{plot the stat.tab file on a map}
\description{
plot the result of on a map
}
\usage{
plotStatMap(df, pop, z.param, margin = 0.1, zlab = z.param, ma = 1, xlim= NULL, ylim=NULL, main=paste(pop, "-population"), track = NULL, cex=1, ...)
}
\arguments{
\item{df}{dataframe of the summary data, ie stat.tab}
\item{pop}{Name should match the one written in the stats.tab created by the Summarize function}
\item{z.param}{parameter in the dataframe for which to plot the heatmap}
\item{zlab}{the label for the values of the heatmap color gradient passed as 'col'}
\item{margin}{margin of longitude and latitude around edges of current position.}
\item{ma}{number of periods to average over z.param}
\item{xlim}{limit of x axis}
\item{ylim}{limit of y axis}
\item{main}{Plot title}
\item{track}{a dataframe of the cruise track (from the sds file) with longitude (long) and latitude (lat) named accordingly.}
\item{cex}{scaling for plotting text and symbols.}
\item{...}{additional arguments to be passed to the plot function}
}

\value{
returns a map of seaflow statistics 
}
\examples{

## load the data
stat.tab <- system.file("extdata","seaflow_cruise","stats.tab", 
					package="flowPhyto")
stats <- read.delim(stat.tab)

## plot the cell concentrations of the picoplankton population
plotStatMap(df=stats ,pop='synecho', z.param='conc')
mtext(line=1, side=4, "cell concentration 10^6 cells / L")
}

