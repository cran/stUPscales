\name{MC.summary}
\alias{MC.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary statistics computation of Monte Carlo simulation}
\description{
A function that computes the summary statistics of a Monte Carlo simulation result.
}
\usage{
MC.summary(p1, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{
The independient variable. A \code{dataframe} object with two columns and number 
of rows equal to the number of rows 
of the Monte Carlo simulated data. The first column, named "time", contains the vector 
of time of the time series in format \code{POSIXct}. The second column 
contains the observations of the time series.
}
  \item{data}{
A \code{matrix} or a \code{dataframe} that contains the results of a Monte Carlo simulation,
with number of rows equal to the number of Monte Carlo realizations and number of columns
equal to the number of oservations i.e. equal to the number of rows of "p1".
}
}
\details{
This function is internally invoked by the \code{MC.analysis} function to compute the summary
statistics of the Monte Carlo simulation under analysis.
}
\value{
A \code{dataframe} with n observations of 15 variables, where n is the number of columns 
of the "data" argument. The 15 variables are time series with the summary statistics of the 
Monte Carlo data: 1) \code{idx}: an index for the dataset 
equal to 1; 2) \code{Mean}: the mean; 3) \code{Sd}: the standard deviation; 
4) \code{Variance}, the variance; 5) \code{q05}: the five percent quantile; 
6) \code{q25}: the 25 percent quantile; 7) \code{q50}: the 50 percent quantile; 
8) \code{q75}: the 75 percent quantile; 9) \code{q95}: the 95 percent quantile; 
10) \code{q995}: the 99.5 percent quantile; 11) \code{q999}: the 99.9 percent quantile; 
12) \code{Max}: the maximum; 13) \code{Sum}: the sum; 14) \code{time}: the time; 
15) \code{p1}: the independient variable.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
J.A. Torres-Matallana
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
library(stUPscales)
library(EmiStatR)

data(P1)
colnames(P1)

new_data <- t(matrix(data = rep(runif(nrow(P1), 10, 100), 5), nrow = nrow(P1), ncol = 5))
new_summary <- MC.summary(p1 = P1, data = new_data)
str(new_summary)
head(new_summary)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Monte Carlo summary statistics}% use one of  RShowDoc("KEYWORDS")
\keyword{MC.summary}% __ONLY ONE__ keyword per line
