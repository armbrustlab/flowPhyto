\name{readSeaflow}
\alias{readSeaflow}
\title{Read a SeaFlow File}
\description{
reads a binary SeaFlow event file and converts into an event dataframe
}
\usage{
readSeaflow(file.path, column.names = EVT.HEADER, column.size = 2, 
		count.only=FALSE,  add.yearday.file=FALSE)
}
\arguments{
\item{file.path}{System path to the binary seaflow event file.}
\item{column.names}{Names of the channels. By default it uses the standard SeaFlow channels described in 'EVT.HEADER' that are
[1] "time" [2] "pulse_width" [3] "D1" [4] "D2" [5] "fsc_small" [6] "fsc_perp" [7] "fsc_big" [8] "pe" [9] "chl_small" [9] "chl_big" }
\item{column.size}{Sizes in bytes of the columns. Set up at 2 by default}
\item{count.only}{Just read the first line of the file (the event/line count). Use to check the conversion of the binary file. FALSE by default}
\item{add.yearday.file}{append the year_day and file number as two final columns in the returned dataframe.}
}

\value{
returns a seaflow dataframe
}
\examples{

opp.file.path <- system.file("extdata","seaflow_cruise","2011_001", "1.evt.opp",
				package="flowPhyto")

opp <- readSeaflow(opp.file.path)
}


