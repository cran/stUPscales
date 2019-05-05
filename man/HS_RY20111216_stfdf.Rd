\name{HS_RY20111216_stfdf}
\alias{HS_RY20111216_stfdf}
\docType{data}
\title{
5-minute DWD precipitation radar imagery non-calibrated in STFDF format
}
\description{
Non-calibrated 5-minute precipitation radar imagery from the German Weather Service
(DWD from the initials in German) at one-kilometer spatial resolution 
over the Haute-Sure catchment of the Grand-Duchy of Luxembourg. 
The data was recorded at Neuheilenbach radar station located in territory of Germany,
which covers the entire territory of the Grand-Duchy of Luxembourg and surroundings.
This sample STFDF (spatio-temporal full data.frame) corresponds to 1-day sample event 
for precipitation recorded on 16th December 2011. 

\emph{Please note that these are un-calibrated radar data.}
}

\usage{data("HS_RY20111216_stfdf")}
\format{
  The format is: \cr
Formal class 'STFDF' [package "spacetime"] with 4 slots \cr
  ..@ data   :'data.frame':	609408 obs. of  1 variable:\cr
  .. ..$ raa01.ry_10000.1112160000.dwd...bin: num [1:609408] 1 2 1 1 1 1 1 1 1 1 ...\cr
  ..@ sp     :Formal class 'SpatialPixels' [package "sp"] with 5 slots\cr
  .. .. ..@ grid       :Formal class 'GridTopology' [package "sp"] with 3 slots\cr
  .. .. .. .. ..@ cellcentre.offset: Named num [1:2] -327 -4347\cr
  .. .. .. .. .. ..- attr(*, "names")= chr [1:2] "s1" "s2"\cr
  .. .. .. .. ..@ cellsize         : num [1:2] 1 1\cr
  .. .. .. .. ..@ cells.dim        : int [1:2] 46 46\cr
  .. .. ..@ grid.index : int [1:2116] 1 2 3 4 5 6 7 8 9 10 ...\cr
  .. .. ..@ coords     : num [1:2116, 1:2] -327 -326 -325 -324 -323 ...\cr
  .. .. .. ..- attr(*, "dimnames")=List of 2\cr
  .. .. .. .. ..$ : NULL\cr
  .. .. .. .. ..$ : chr [1:2] "s1" "s2"\cr
  .. .. ..@ bbox       : num [1:2, 1:2] -327 -4348 -281 -4302\cr
  .. .. .. ..- attr(*, "dimnames")=List of 2\cr
  .. .. .. .. ..$ : chr [1:2] "s1" "s2"\cr
  .. .. .. .. ..$ : chr [1:2] "min" "max"\cr
  .. .. ..@ proj4string:Formal class 'CRS' [package "sp"] with 1 slot\cr
  .. .. .. .. ..@ projargs: chr "+proj=stere +lat_0=90 +lat_ts=90 +lon_0=10 +k=0.93301270189 +x_0=0 +y_0=0 +a=6370040 +b=6370040 +units=km +no_defs"\cr
  ..@ time   :An 'xts' object on 2011-12-16/2011-12-16 23:55:00 containing:\cr
  Data: int [1:288, 1] 1 2 3 4 5 6 7 8 9 10 ...\cr
 - attr(*, "dimnames")=List of 2\cr
  ..$ : NULL\cr
  ..$ : chr "timeIndex"\cr
  Indexed by objects of class: [POSIXct,POSIXt] TZ: \cr
  xts Attributes:  \cr
 NULL\cr
  ..@ endTime: POSIXct[1:288], format: "2011-12-16 00:00:00" "2011-12-16 00:05:00" "2011-12-16 00:10:00" "2011-12-16 00:15:00" ...\cr
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
library(spacetime)

data(HS_RY20111216_stfdf)

sample.idx <- seq.default(from = 1, to = 25, by = 1)
sample     <- HS_RY20111216_stfdf[, sample.idx]

stplot(sample,
       scales=list(draw=TRUE),
       key.space="right", colorkey=TRUE,
       main="5-minute DWD sample precipitation radar imagery non-calibrated in STFDF format",
       cex=.74, par.strip.text=list(cex=.74))
}

\keyword{datasets}
