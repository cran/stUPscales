\name{Agg.t}
\alias{Agg.t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Temporal aggregation of environmental variables
}
\description{
  Function for temporal aggregation of environmental variables. \code{Agg} is a wrapper function of 
  \code{aggregate} from \code{stats} package.
}
\usage{
  Agg.t(data, nameData, delta, func, namePlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    A \code{data.frame} that contains the time series of the 
    environmental variable to be aggregated, e.g. precipitation. This data.frame should have 
    at two columns: the first one, Time [y-m-d h:m:s]; the second one, a numeric value equal to the 
    magnitude of the environmental variable. If the environmental variable is different than precipitation,
    then the column name of the values can be named as the name of the variable itself.}
  \item{nameData}{
    A \code{character} string that defines the name of the environmental variable to be aggregated.}
  \item{delta}{
    A \code{numeric} value that specifies the level of aggregation required in minutes.}
  \item{func}{
    The name of the function of aggregation e.g. \code{mean}, \code{sum}.}
  \item{namePlot}{
    A \code{character} string that defines the title of the plot generated.}
}
%\details{
  %%  ~~ If necessary, more details than the description above ~~
    %}
\value{
  A \code{data.frame} with two columns:
    \item{time}{the date-time time series of the aggregated variable}
    \item{value}{time series with the magnitude of the aggregated variable.}
}
%\references{
  %% ~put references to the literature/web site here ~
    %}
\author{
  J.A. Torres-Matallana
}
%\note{
  %  
  %}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
%  \seealso{
%    See also \code{\link{setup}}, \code{\link{MC.setup}}, \code{\link{MC.sim}}.
%  }

\examples{

## temporal aggregation
library(EmiStatR)
data(P1)
colnames(P1) <- c("time", "P1")
head(P1)

library(stUPscales)
P1.agg <- Agg.t(data = P1, nameData = "P1", delta = 120 , func = sum, 
              namePlot = "Temporal aggregation of precipitation P1")

head(P1.agg)
tail(P1.agg)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Temporal aggregation}
\keyword{Agg.t}
\keyword{Aggregation}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line