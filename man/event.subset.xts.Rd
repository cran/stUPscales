\name{Lux_precipitation}
\alias{Lux_precipitation}
\alias{event.subset.xts}
\docType{data}
\title{
Sample precipitation time series in the Grand-Duchy of Luxembourg
}
\description{
A 10-hour sample event for precipitation time series measured in 25 rain gauge stations
distributed over the territory of the Grand-Duchy of Luxembourg.
}
\usage{data("Lux_precipitation")}
\format{
  The format is: \cr
An 'xts' object on 2011-12-16/2011-12-16 10:00:00 containing:\cr
  Data: num [1:61, 1:25] 0 0 0 0 0.1 0.1 0.1 0.1 0.1 0.2 ...\cr
 - attr(*, "dimnames")=List of 2\cr
  ..$ : NULL\cr
  ..$ : chr [1:25] "Dahl" "Echternach" "Esch-Sure" "Eschdorf" ...\cr
  Indexed by objects of class: [POSIXct,POSIXt] TZ: \cr
  xts Attributes:  \cr
 NULL\cr
}
%\details{
%%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
\url{http://agrimeteo.lu/}
}
%\references{
%%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(Lux_precipitation)

library(xts)
head(event.subset.xts)
tail(event.subset.xts)

plot(event.subset.xts)
}
\keyword{datasets}
