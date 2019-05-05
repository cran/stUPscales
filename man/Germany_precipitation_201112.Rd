\name{Germany_precipitation_201112}
\alias{Germany_precipitation_201112}
\docType{data}
\title{
Sample precipitation time series in Germany
}
\description{
A 1-minute sample event for precipitation time series measured in 37 rain gauge stations
distributed over the territory of Germany close to the frontier to the Grand-Duchy of Luxembourg.
}
\usage{data("Germany_precipitation_201112")}
\format{
  The format is:\cr
An 'xts' object on 2011-12-01/2011-12-31 23:59:00 containing:\cr
  Data: num [1:44640, 1:37] NA NA NA NA NA NA NA NA NA NA ...\cr
 - attr(*, "dimnames")=List of 2\cr
  ..$ : NULL\cr
  ..$ : chr [1:37] "Station.1327" "Station.1964" "Station.200" "Station.2170" ...\cr
  Indexed by objects of class: [POSIXct,POSIXt] TZ: \cr
  xts Attributes:  \cr
 NULL\cr
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
\url{https://www.dwd.de/}
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
library(stUPscales)

data(Germany_precipitation_201112)

summary(Germany_precipitation_201112)

library(zoo)
par(mfrow = c(5, 1))
for(i in 1:ncol(Germany_precipitation_201112)){
  plot(index(Germany_precipitation_201112), Germany_precipitation_201112[,i], 
       typ="l", col="blue", xlab="Time", ylab='Precip. [mm]', 
       main=colnames(Germany_precipitation_201112[,i]))
}

}
\keyword{datasets}
