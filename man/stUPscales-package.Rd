\name{stUPscales-package}
\alias{stUPscales-package}
\alias{stUPscales}
\docType{package}
\title{Spatio-Temporal Uncertainty Propagation Across Multiple Scales}
\description{Integrated environmental modelling requires coupling sub-models at different spatial and temporal scales, thus accounting for change of support procedures (aggregation and disaggregation). We contribute to state-of-the-art open source tools that support uncertainty propagation analysis in temporal and spatio-temporal domains. We implement the tool for uncertainty propagation in environmental modelling, with examples in the urban water domain. The functionalities of the class setup and the methods and functions MC.setup, MC.sim, MC.analysis, MC.analysis_generic and Agg.t are contained, which are used for setting up, running and analysing Monte Carlo uncertainty propagation simulations, and for spatio-temporal aggregation. We also implement functionalities to model and predict variables that vary in space and time. stUPscales takes uncertainty characterisation and propagation a step further by including temporal and spatio-temporal auto- and cross-correlation, resulting in more realistic (spatio-)temporal series of environmental variables. Due to its modularity, the package allows the implementation of additional methods and functions for spatio-temporal disaggregation of model inputs and outputs, when linking models across multiple space-time scales.     
}
\details{

The DESCRIPTION file:
\tabular{ll}{
Package: \tab stUPscales  \cr
Type:    \tab Package     \cr
Version: \tab 1.0.3.4     \cr
Date:    \tab 2019-05-04  \cr
License: \tab GPL (>= 3)  \cr
Depends: \tab R (>= 2.10), methods, stats, graphics, grDevices, utils, mAr, lmom \cr
Imports: \tab parallel, doParallel, foreach, lattice, msm, ggplot2, moments, hydroGOF, zoo, data.table, xts, EmiStatR \cr
Suggests: \tab sp, spacetime \cr
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
