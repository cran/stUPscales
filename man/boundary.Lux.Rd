\name{Lux_boundary}
\alias{Lux_boundary}
\alias{boundary.Lux}
\docType{data}
\title{
A shapefile for the boundary of the Grand-Duchy of Luxembourg
}
\description{
A shapefile for the country boundary of the Grand-Duchy of Luxembourg
}
\usage{data("Lux_boundary")}
\format{
  The format is:\cr
Formal class 'SpatialPolygonsDataFrame' [package "sp"] with 5 slots\cr
  ..@ data       :'data.frame':	1 obs. of  3 variables:\cr
  .. ..$ cat: int 1\cr
  .. ..$ X_ : Factor w/ 1 level "?": 1\cr
  .. ..$ X_1: Factor w/ 1 level "?": 1\cr
  ..@ polygons   :List of 1\cr
  .. ..$ :Formal class 'Polygons' [package "sp"] with 5 slots\cr
  .. .. .. ..@ Polygons :List of 1\cr
  .. .. .. .. ..$ :Formal class 'Polygon' [package "sp"] with 5 slots\cr
  .. .. .. .. .. .. ..@ labpt  : num [1:2] 74692 93669\cr
  .. .. .. .. .. .. ..@ area   : num 2.6e+09\cr
 .. .. .. .. .. .. ..@ hole   : logi FALSE\cr
  .. .. .. .. .. .. ..@ ringDir: int 1\cr
  .. .. .. .. .. .. ..@ coords : num [1:1141, 1:2] 82274 82437 82449 82765 82782 ...\cr
  .. .. .. ..@ plotOrder: int 1\cr
  .. .. .. ..@ labpt    : num [1:2] 74692 93669\cr
  .. .. .. ..@ ID       : chr "0"\cr
  .. .. .. ..@ area     : num 2.6e+09\cr
  ..@ plotOrder  : int 1\cr
  ..@ bbox       : num [1:2, 1:2] 49034 57132 106245 138879\cr
  .. ..- attr(*, "dimnames")=List of 2\cr
  .. .. ..$ : chr [1:2] "x" "y"\cr
  .. .. ..$ : chr [1:2] "min" "max"\cr
  ..@ proj4string:Formal class 'CRS' [package "sp"] with 1 slot\cr
  .. .. ..@ projargs: chr "+proj=tmerc +lat_0=49.83333333333334 +lon_0=6.166666666666667 +k=1 +x_0=80000 +y_0\%\%=100000 +ellps=intl +units=m +no_defs"\cr
}
%\details{
%%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
library(stUPscales)

data(Lux_boundary)

str(boundary.Lux)
}
\keyword{datasets}
