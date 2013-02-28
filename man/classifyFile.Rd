\name{classifyFile}
\alias{classifyFile}
\title{Cluster the different Phytoplankton Populations}
\description{
classify the different cell populations from an OPP or FCS file according to a  pre-defined parameters of population definition and output a consensus vector.

For each group of OPP or FCS file, the function outputs a single vector file (consensus.vct) that contains the population identification of each single cell.
The function is run in single OPP or FCS file increments to provide multiple vector files over each OPP or FCS file and strengthen the clustering analysis.

}
\usage{
 classifyFile(opp.path, concat.ct=3, output.path=getCruisePath(opp.path), 
 		def.path=paste(getCruisePath(opp.path),'pop.def.tab',sep=''), varnames= CHANNEL.CLMNS.SM, numc=0, noise=0, auto.correction=FALSE)
}
\arguments{
\item{opp.path}{Path to OPP file.}
\item{concat.ct}{Number of OPP files to concatenate at a time.}
\item{output.path}{Path to the directory where you wish to output data.}
\item{varnames}{A character vector specifying the variables (columns) to be included in clustering when choosing flowMeans.}
\item{numc}{Number of clusters when choosing flowMeans. If set to 0 (default) the value matches the number of populations defined in pop.def table . If set to NA, the optimal number of clusters will be estimated automatically.}
\item{noise}{Set up the noise threshold for phytoplankton cells. Only cells with chlorophyll value higher than the noise will be clustered}
\item{auto.correction}{Correct for potential mislabeling of phytoplankton populations defined in pop.def.tab}
\item{def.path}{Path to the file that defines how to gate & cluster the events into populations.}
}
\value{
a consensus vector file composed of one single column indicating population assignment of each event from the OPP or FCS file.
}

\examples{
example.cruise.name <- 'seaflow_cruise'
temp.out.dir <- '.'

seaflow.path <- system.file("extdata", example.cruise.name, package="flowPhyto")
file.copy(from=seaflow.path, to=temp.out.dir, recursive=TRUE)

classifyFile(paste(temp.out.dir,'/',example.cruise.name,'/2011_001/1.evt.opp',sep=''))
unlink(example.cruise.name, recursive=TRUE)

}


