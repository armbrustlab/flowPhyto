\name{pipeline}
\alias{pipeline}
\title{Run the SeaFlow Pipeline}
\description{
run the pipeline
}
\usage{
pipeline(cruise.name='', repo=REPO.PATH, range=NULL, steps=1:4, pct=.97,  
	clust.concat.ct=3, resample.size=300, resamp.concat.max=10,
	 filter.width=1.5, filter.notch=1,  filter.edge=1, filter.slope=NA,
	 classify.func =2, classify.varnames=CHANNEL.CLMNS.SM, classify.numc=0, classify.noise=0,
	 map.margin=2, 
	concat.sds=!is.na(match(1,steps)), load.to.db=FALSE,  preplot=FALSE, cleanup=TRUE,
	input.path=paste(repo, '/', cruise.name, sep=''),
	output.path=input.path,  log.dir=output.path,
	def.path=paste(input.path,'/', 'pop.def.tab',sep=''), parallel=TRUE, submit.cmd='qsub')
}
\arguments{
\item{cruise.name}{Simplified cruise name (same name as the subdirectory in the seaflow data dir).}
\item{steps}{Which steps of the pipeline to run. step 1 is filter, step 2 is classify, step 3 is census and consensus, step 4 is summarize. 1:2 will do step 1 to 2, etc.}
\item{pct}{percentage completion (number of indicator files created vs input files) each job step should go to.}
\item{clust.concat.ct}{Number of event file to concatenate at a time during the clustering/classification step.}
\item{map.margin}{Margin in latitude/longitude around the map plots.}
\item{resample.size}{Minimum number of events in a population.}
\item{resamp.concat.max}{Maximum number of allowable event files to concatenate to generate statistics from.}
\item{filter.notch}{the location of the x=y (by default) point to create the notch in the gated filter}
\item{filter.width}{the margin of error for particle alignment determination in the filter step.}
\item{filter.slope}{correction factor for the stream alignment. When stream is not properly aligned, aligned particles do not scatter light equally on D1 and D2 and therefore do not lie onto the 1:1 line on the scatter plot of D2 vs D1. By default, the value of the slope is calculated as the ratio of D2/D1}
\item{filter.edge}{location of the boundary layer between water/air. Particles located at the boundary layer scatter light that can be detected by the position detectors for the filter step.}
\item{classify.func}{Choose the clustering method, either flowClust (func = 1) or flowMeans (func = 2, by default) function}
\item{classify.varnames}{A character vector specifying the variables (columns) to be included in clustering when choosing flowMeans.}
\item{classify.numc}{Number of clusters when choosing flowMeans. If set to 0 (default) the value matches the number of populations defined in pop.def table . If set to NA, the optimal number of clusters will be estimated automatically.}
\item{classify.noise}{Set up the noise threshold for phytoplankton cells. Only cells with chlorophyll value higher than the noise will be clustered}
\item{concat.sds}{Determines if the sds files in the individual julian day directories should be concatenated together into sds.tab}
\item{load.to.db}{Load the sds and stat files to the database.}
\item{preplot}{Preplot the level 2 analysis plots to 'output.path'.}
\item{cleanup}{Cleanup the submission and (non error reporting) R CMD BATCH log files.}
\item{input.path}{Path to the directory with input data (raw evt or opp files.}
\item{output.path}{Path to the directory where you wish to output data.}
\item{log.dir}{Path to the directory where log file will be written.}
\item{def.path}{Path to the file that defines how to gate & cluster the events into populations.}
\item{parallel}{Boolean indicating if the job should be run in parallel using qsub (vs in serial)}
\item{repo}{Full path to your SeaFlow repository}
\item{range}{A named, two-integer vector specifying the start and end (inclusive) range for subsetting the input files used in each analysis step (with the exception of summarize). Values should be a (evt/opp) file numbers and names should be strings corresponding to the year_julianday directory names.  The nv() function is useful for creating this vector.}
\item{submit.cmd}{the command used to deploy an R CMD BATCH system call to a cluster.  Must be used in conjunction with parallel=TRUE.}
}

\examples{

example.cruise.name <- 'seaflow_cruise'
temp.out.dir <- '.' #path.expand('~')
output.path <- paste(temp.out.dir,'/',example.cruise.name,sep='')
seaflow.path <- system.file("extdata", example.cruise.name, package="flowPhyto")

file.copy(from=seaflow.path, to=temp.out.dir, recursive=TRUE)

pipeline(repo= temp.out.dir, cruise.name='seaflow_cruise', steps=4, parallel=FALSE) 
unlink(example.cruise.name, recursive=TRUE)


}


