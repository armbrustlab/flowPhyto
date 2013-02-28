\name{consensus}
\alias{consensus}
\title{create a consensus vector from several classification vectors}
\description{
Create a consensus from several classification vectors 
}
\usage{
consensus(mtrx, threshold=.6)
}
\arguments{
\item{mtrx}{a matrix of n concatenated classification vectors.}
\item{threshold}{the minimum percentage of identical classification calls required for an unambiguous consensus call.}
}
\value{
a vector the same length as the number of rows in the input matrix with a population classification call for each element
}

\examples{
vct.paths <- sapply(c(1,439,440), function(i) 
		system.file("extdata","seaflow_cruise","2011_001", 
			paste("1.evt.opp.",i,'-class.vct',sep=''), 
				package="flowPhyto"))
mat <- do.call(cbind,lapply(vct.paths, read.delim))
v <- consensus(mtrx=mat)
table(v$pop)
aggregate(v$support,list(v$pop), mean)
}


