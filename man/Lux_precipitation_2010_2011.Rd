\name{Lux_precipitation_2010_2011}
\alias{Lux_precipitation_2010_2011}
\docType{data}
\title{
Sample precipitation time series in the Grand-Duchy of Luxembourg (2-year period)
}
\description{
A 2-year period sample for precipitation time series measured at 10-minute time step
in 25 rain gauge stations distributed over the territory of the Grand-Duchy of Luxembourg.
}
\usage{data("Lux_precipitation_2010_2011")}
\format{
  The format is:\cr
An 'xts' object on 2010-01-01/2011-12-31 23:50:00 containing:\cr
  Data: num [1:105120, 1:25] 0 0 0 0 0 0 0 0 0 0 ...\cr
 - attr(*, "dimnames")=List of 2\cr
  ..$ : NULL\cr
  ..$ : chr [1:25] "Dahl" "Echternach" "Esch-Sure" "Eschdorf" ...\cr
  Indexed by objects of class: [POSIXct,POSIXt] TZ: \cr
  xts Attributes:  \cr
 NULL\cr
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
\url{https://www.agrimeteo.lu/}
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
library(stUPscales)
data(Lux_precipitation_2010_2011)

library(xts)
head(Lux_precipitation_2010_2011)
tail(Lux_precipitation_2010_2011)

## uncomment to run (plot)
# plot(Lux_precipitation_2010_2011)
}
\keyword{datasets}
