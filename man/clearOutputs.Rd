\name{clearOutputs}
\alias{clearOutputs}
\title{Clear the output files from particular flowPhyto pipeline steps.}
\description{
Each step of the flowPhyto pipeline generates many files which are used as indicators of completion. This function helps to clear those files away to allow for proper waiting between steps of a pipeline rerun.
}
\usage{
clearOutputs(cruise.path='.', steps=1:4)
}
\arguments{
\item{cruise.path}{Path to the cruise directory.}
\item{steps}{Steps for which to clear the outputs.}
}
\value{
none
}
\examples{

example.cruise.name <- 'seaflow_cruise'
temp.out.dir <- '.'

seaflow.path <- system.file("extdata", example.cruise.name, package="flowPhyto")
file.copy(from=seaflow.path, to=temp.out.dir, recursive=TRUE)

clearOutputs(paste(temp.out.dir,'/',example.cruise.name,sep=''), steps=3)
unlink(example.cruise.name, recursive=TRUE)

}

