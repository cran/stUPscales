\name{inputObs-class}
\Rdversion{1.1}
\docType{class}
\alias{inputObs-class}
\alias{inputObs}
%\alias{Validation_Quantity_Agg,input,inputObs-method}
%\alias{Validation_Quantity,input,inputObs-method}

\title{Class \code{"inputObs"}}
\description{
The class provides a container for inputs required to invoke 
\code{Validation_Quantity} and \code{Validation_Quantity_Agg} methods.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{inputObs(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{id}:}{Object of class \code{"numeric"} to define an unique index for the object.}
    \item{\code{plot}:}{Object of class \code{"numeric"}. One of 0 (no plots are cretated) or 
    1 (to create plots).}
    \item{\code{delta}:}{Object of class \code{"list"} to define the time step in minutes for 
    temporal aggregation required e.g. list(P1 = 10, wlt_obs = 10, vol_sim = 10)} 
    for defining the time steps of 10 minutes for the three variables P1, wlt_obs, vol_sim.
    \item{\code{observations}:}{Object of class \code{"list"} to define the observed time series.}
    \item{\code{lev2vol}:}{Object of class \code{"list"} to define the curve for the 
    relationship level to volume.}
    \item{\code{namePlot}:}{Object of class \code{"character"} to define the name of the plot to create.}
    \item{\code{legendPosition}:}{Object of class \code{"list"} with three character objects, 
    which define the posistion of the legend for the top, second and bottom insets of the plot. }
    \item{\code{var}:}{Object of class \code{"character"} to define the name of the variable 
    from which the time series simulated are defined. }
  }
}
\section{Methods}{
  \describe{
    \item{Validation_Quantity_Agg}{\code{signature(x = "input", y = "inputObs")}: ... }
    \item{Validation_Quantity}{\code{signature(x = "input", y = "inputObs")}: ... }
	 }
}
%\references{
%%%  ~~put references to the literature/web site here~~
%}
\author{
J.A. Torres-Matallana
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("inputObs")

inputObs()
}
\keyword{classes}

%\usage{
%inputObs(...)
%}

%\arguments{
%  \item{\dots}{
%%%     ~~Describe \code{\dots} here~~
%}
%}
