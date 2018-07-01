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
  The format is:
Formal class 'SpatialPolygonsDataFrame' [package "sp"] with 5 slots
  ..@ data       :'data.frame':	1 obs. of  3 variables:
  .. ..$ cat: int 1
  .. ..$ X_ : Factor w/ 1 level "?": 1
  .. ..$ X_1: Factor w/ 1 level "?": 1
  ..@ polygons   :List of 1
  .. ..$ :Formal class 'Polygons' [package "sp"] with 5 slots
  .. .. .. ..@ Polygons :List of 1
  .. .. .. .. ..$ :Formal class 'Polygon' [package "sp"] with 5 slots
  .. .. .. .. .. .. ..@ labpt  : num [1:2] 74692 93669
  .. .. .. .. .. .. ..@ area   : num 2.6e+09
 .. .. .. .. .. .. ..@ hole   : logi FALSE
  .. .. .. .. .. .. ..@ ringDir: int 1
  .. .. .. .. .. .. ..@ coords : num [1:1141, 1:2] 82274 82437 82449 82765 82782 ...
  .. .. .. ..@ plotOrder: int 1
  .. .. .. ..@ labpt    : num [1:2] 74692 93669
  .. .. .. ..@ ID       : chr "0"
  .. .. .. ..@ area     : num 2.6e+09
  ..@ plotOrder  : int 1
  ..@ bbox       : num [1:2, 1:2] 49034 57132 106245 138879
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:2] "x" "y"
  .. .. ..$ : chr [1:2] "min" "max"
  ..@ proj4string:Formal class 'CRS' [package "sp"] with 1 slot
  .. .. ..@ projargs: chr "+proj=tmerc +lat_0=49.83333333333334 +lon_0=6.166666666666667 +k=1 +x_0=80000 +y_0\%\%=100000 +ellps=intl +units=m +no_defs"
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
data(Lux_boundary)

str(boundary.Lux)
}
\keyword{datasets}
