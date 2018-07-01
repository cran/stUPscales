\name{Lux_stations}
\alias{Lux_stations}
\alias{stations}
\docType{data}
\title{
A shapefile for the location of 25 rain gauges
}
\description{
A shapefile for the location of 25 rain gauges distributed over the territory of the Grand-Duchy of Luxembourg.
This 25 stations are the same related in the "event.subset.xts" dataset.

}
\usage{data("Lux_stations")}
\format{
  The format is:
Formal class 'SpatialPointsDataFrame' [package "sp"] with 5 slots
  ..@ data       :'data.frame':	25 obs. of  8 variables:
  .. ..$ id        : Factor w/ 25 levels "1","11","12",..: 4 5 6 7 8 9 10 11 25 1 ...
  .. ..$ name      : Factor w/ 25 levels "Arsdorf","Christnach",..: 4 5 7 6 8 9 10 11 22 1 ...
  .. ..$ north_lure: Factor w/ 24 levels "101950","102913",..: 6 NA 4 3 1 5 22 20 21 2 ...
  .. ..$ east_luref: Factor w/ 25 levels "56584","56990",..: 10 25 5 7 16 19 20 24 17 1 ...
  .. ..$ elev_luref: Factor w/ 25 levels "190","202","207",..: 22 5 16 25 2 11 14 1 6 20 ...
  .. ..$ station_ty: Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ management: Factor w/ 1 level "ASTA": 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ telemetry : Factor w/ 0 levels: NA NA NA NA NA NA NA NA NA NA ...
  ..@ coords.nrs : num(0) 
  ..@ coords     : num [1:25, 1:2] 66562 99810 62258 63363 74929 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:2] "coords.x1" "coords.x2"
  ..@ bbox       : num [1:2, 1:2] 56584 64215 99810 132012
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:2] "coords.x1" "coords.x2"
  .. .. ..$ : chr [1:2] "min" "max"
  ..@ proj4string:Formal class 'CRS' [package "sp"] with 1 slot
  .. .. ..@ projargs: chr "+proj=tmerc +lat_0=49.83333333333334 +lon_0=6.166666666666667 +k=1 +x_0=80000 +y_0=100000 +ellps=intl +units=m +no_defs"
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
data(Lux_stations)

str(stations)

}
\keyword{datasets}
