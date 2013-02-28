\name{classify}
\alias{classify}
\title{Cluster the different Phytoplankton Populations}
\description{
Classify the different cell populations from an OPP or FCS dataframe according to a  pre-defined parameters of population definition

Because the characteristics of each phytoplankton populations varied according to environmental conditions and instrument settings, a customizable table of pre-defined parameters (pop.def) is used to help in gating the different phytoplankton populations. 
The rows of the pop.def table represent the names of the different populations.
The columns of the pop.def table represent the parameters used for gating and clustering the different populations.
The function uses these pre-defined parameters and inputs a single OPP or FCS file to cluster cell populations using flowMeans package 
}
\usage{
classify(x, pop.def=POP.DEF, varnames = CHANNEL.CLMNS.SM, numc=0, noise=0,  auto.correction=FALSE, plot.cluster=FALSE, plot.assignment=FALSE, try.context='local',...)
}
\arguments{
\item{x}{an OPP or FCS dataframe.}
\item{pop.def}{pop.def table that defines how to gate & cluster the events into populations.}
\item{varnames}{A character vector specifying the variables (columns) to be included in clustering when choosing flowMeans.}
\item{numc}{Number of clusters when choosing flowMeans. If set to 0 (default) the value matches the number of populations defined in pop.def table . If set to NA, the optimal number of clusters will be estimated automatically.}
\item{noise}{Set up the noise threshold for phytoplankton cells. Only cells with chlorophyll value higher than the noise will be clustered}
\item{auto.correction}{Correct for potential mislabeling of phytoplankton populations defined in pop.def.tab}
\item{plot.cluster}{Plot the output of clustering when choosing flowMeans}
\item{plot.assignment}{Plot the output of Matching cluster number with cell population defined in pop.def.tab when choosing flowMeans}
\item{try.context}{Default value set up to 'local'}
\item{...}{additional arguments to be passed to the plot function}

}
\value{
an OPP or FCS dataframe like the input x but with an additional column 'pop' indicating population assignment
}
\examples{

## reading a standard SeaFlow file
opp.path <- system.file("extdata","seaflow_cruise","2011_001", "1.evt.opp", 
				package="flowPhyto")
pop.def.path <- system.file("extdata","seaflow_cruise","pop.def.tab", 
				package="flowPhyto")

opp <- readSeaflow(opp.path)
def <- readPopDef(pop.def.path)
pop <- classify(x=opp, pop.def= def)

table(pop$pop)

}


