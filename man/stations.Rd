\name{Lux_stations}
\alias{Lux_stations}
\alias{stations}
\docType{data}
\title{
A SpatialPointsDataFrame with the location of 25 rain gauges in Luxembourg
}
\description{
A SpatialPointsDataFrame with the location of 25 rain gauges distributed over the territory of the Grand-Duchy of Luxembourg.
These 25 stations are the same related to the "event.subset.xts" dataset.

}
\usage{data("Lux_stations")}
\format{
  The format is:\cr
Formal class 'SpatialPointsDataFrame' [package "sp"] with 5 slots\cr
  ..@ data       :'data.frame':	25 obs. of  8 variables:\cr
  .. ..$ id        : Factor w/ 25 levels "1","11","12",..: 4 5 6 7 8 9 10 11 25 1 ...\cr
  .. ..$ name      : Factor w/ 25 levels "Arsdorf","Christnach",..: 4 5 7 6 8 9 10 11 22 1 ...\cr
  .. ..$ north_lure: Factor w/ 24 levels "101950","102913",..: 6 NA 4 3 1 5 22 20 21 2 ...\cr
  .. ..$ east_luref: Factor w/ 25 levels "56584","56990",..: 10 25 5 7 16 19 20 24 17 1 ...\cr
  .. ..$ elev_luref: Factor w/ 25 levels "190","202","207",..: 22 5 16 25 2 11 14 1 6 20 ...\cr
  .. ..$ station_ty: Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ management: Factor w/ 1 level "ASTA": 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ telemetry : Factor w/ 0 levels: NA NA NA NA NA NA NA NA NA NA ...\cr
  ..@ coords.nrs : num(0) \cr
  ..@ coords     : num [1:25, 1:2] 66562 99810 62258 63363 74929 ...\cr
  .. ..- attr(*, "dimnames")=List of 2\cr
  .. .. ..$ : NULL\cr
  .. .. ..$ : chr [1:2] "coords.x1" "coords.x2"\cr
  ..@ bbox       : num [1:2, 1:2] 56584 64215 99810 132012\cr
  .. ..- attr(*, "dimnames")=List of 2\cr
  .. .. ..$ : chr [1:2] "coords.x1" "coords.x2"\cr
  .. .. ..$ : chr [1:2] "min" "max"\cr
  ..@ proj4string:Formal class 'CRS' [package "sp"] with 1 slot\cr
  .. .. ..@ projargs: chr "+proj=tmerc +lat_0=49.83333333333334 +lon_0=6.166666666666667 +k=1 +x_0=80000 +y_0=100000 +ellps=intl +units=m +no_defs"\cr
}
%\details{
%%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
\url{https://www.agrimeteo.lu/}
}
%\references{
%%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
library(stUPscales)
data(Lux_stations)

str(stations)

library(sp)
plot(stations)

}
\keyword{datasets}
