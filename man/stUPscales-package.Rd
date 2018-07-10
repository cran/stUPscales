\name{stUPscales-package}
\alias{stUPscales-package}
\alias{stUPscales}
\docType{package}
\title{Spatio-Temporal Uncertainty Propagation Across Multiple Scales}
\description{Provides several R functions for 
    temporal aggregation of environmental variables used in e.g. Urban Drainage Models (UDMs), 
    as precipitation and pollutants.
    Also, it provides methods and functions for uncertainty propagation via Monte Carlo simulation. 
    This package, moreover, provides specific analysis functions for urban drainage 
    system simulation to evaluate water quantity and quality in combined sewer overflows (CSOs).     
}
\details{

The DESCRIPTION file:
\tabular{ll}{
Package: \tab stUPscales  \cr
Type:    \tab Package   \cr
Version: \tab 1.0.3.2       \cr
Date:    \tab 2018-07-04\cr
License: \tab GPL (>= 3)\cr
Depends: \tab R (>= 2.10), methods, stats, graphics, grDevices, utils, mAr, lmom \cr
Imports: \tab parallel, doParallel, foreach, lattice, msm, ggplot2, moments, hydroGOF, zoo, data.table, xts, EmiStatR \cr
}
}
\author{
    J.A. Torres-Matallana [aut, cre];
    U. Leopold [ctb];
    G.B.M. Heuvelink [ctb].

Maintainer: J.A. Torres-Matallana.
}

%\references{
%~~ Literature or other references for background information ~~
%}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%~~ simple examples of the most important functions ~~
%}
