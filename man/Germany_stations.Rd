\name{Germany_stations}
\alias{Germany_stations}
\docType{data}
\title{
A SpatialPointsDataFrame with the location of 37 rain gauges in Germany
}
\description{
A SpatialPointsDataFrame with the location of 37 rain gauges distributed over the territory of Germany close
to frontier with the Grand-Duchy of Luxembourg.
These 37 stations are the same related to the "Germany_precipitation_201112" dataset.
}
\usage{data("Germany_stations")}
\format{
  The format is:\cr
Formal class 'SpatialPointsDataFrame' [package "sp"] with 5 slots\cr
  ..@ data       :'data.frame':	37 obs. of  9 variables:\cr
  .. ..$ Stations_id  : int [1:37] 200 450 460 523 603 723 902 942 953 1327 ...\cr
  .. ..$ von_datum    : int [1:37] 20020924 20050920 19930930 20020807 20071024 20020717 20060618 20020925 19970730 20040707 ...\cr
  .. ..$ bis_datum    : int [1:37] 20180820 20121204 20180820 20180604 20180820 20180820 20180820 20180820 20180820 20180820 ...\cr
  .. ..$ Stationshoehe: int [1:37] 517 120 363 359 159 290 573 308 481 147 ...\cr
  .. ..$ geoBreite    : num [1:37] 50.1 49.9 49.3 50 50.7 ...\cr
  .. ..$ geoLaenge    : num [1:37] 6.32 7.07 6.69 6.53 7.19 ...\cr
  .. ..$ Stationsname : Factor w/ 1109 levels "Aachen","Aachen-Orsbach",..: 535 90 92 104 118 456 172 179 180 1046 ...\cr
  .. ..$ Bundesland   : Factor w/ 16 levels "Baden-Wuerttemberg",..: 11 11 12 11 10 11 10 11 11 10 ...\cr
  .. ..$ d            : logi [1:37] NA NA NA NA NA NA ...\cr
  ..@ coords.nrs : num(0) \cr
  ..@ coords     : num [1:37, 1:2] 90590 144511 117752 106135 152358 ...\cr
  .. ..- attr(*, "dimnames")=List of 2\cr
  .. .. ..$ : NULL\cr
  .. .. ..$ : chr [1:2] "coords.x1" "coords.x2"\cr
  ..@ bbox       : num [1:2, 1:2] 82780 31291 156554 200669\cr
  .. ..- attr(*, "dimnames")=List of 2\cr
  .. .. ..$ : chr [1:2] "coords.x1" "coords.x2"\cr
  .. .. ..$ : chr [1:2] "min" "max"\cr
  ..@ proj4string:Formal class 'CRS' [package "sp"] with 1 slot\cr
  .. .. ..@ projargs: chr "+init=epsg:2169 +proj=tmerc +lat_0=49.83333333333334 +lon_0=6.166666666666667 +k=1 +x_0=80000 +y_0=100000 +ellp"| __truncated__\cr
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
library(sp)

data(Germany_stations)

str(Germany_stations)

data(Lux_boundary)
plot(Germany_stations)
plot(boundary.Lux, add=TRUE) # Luxembourg boundary
}
\keyword{datasets}
