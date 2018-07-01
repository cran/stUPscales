\name{MC.summary.agg}
\alias{MC.summary.agg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary statistics computation of aggregated Monte Carlo simulation
}
\description{
A function that computes the summary statistics of aggregated Monte Carlo simulation result.
}
\usage{
MC.summary.agg(summ, det, delta, func.agg, func.agg.p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summ}{
A \code{dataframe} with n observations of 15 variables, where n is the number observations 
or time steps of the data. The 15 variables are time series with the summary statistics of the 
Monte Carlo data. This dataframe is in the format as is described in the 
\code{MC.summary} function value.
}
  \item{det}{
A \code{dataframe} that contains the deterministic simulation.}
  \item{delta}{
A \code{numeric} value that represents the level of aggregation (required time stemp) in minutes.
}
  \item{func.agg}{
The aggregation  function to be applied to the \code{summ} dataframe.
}
  \item{func.agg.p}{
The aggregation function to be applied to the independient variable \code{p1} from \code{summ} dataframe.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A \code{dataframe} containing the \code{summ} data aggregated to the level defined by \code{delta}
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

\seealso{
See Also as \code{\link{MC.summary}}
}
\examples{
library(stUPscales)
library(EmiStatR)

data(P1)
colnames(P1)

new_data <- t(matrix(data = rep(runif(nrow(P1), 10, 100), 5), nrow = nrow(P1), ncol = 5))
new_summary <- MC.summary(p1 = P1, data = new_data)
str(new_summary)
head(new_summary)

# deterministic simulation
det <- rnorm(nrow(P1), 45, .15)

# level of aggregation
delta <- 60*2 # 2 hours

new_summary_agg <- MC.summary.agg(summ = new_summary, det, delta, func.agg = mean, func.agg.p = sum)
str(new_summary_agg)
head(new_summary_agg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MC.summary.agg}% use one of  RShowDoc("KEYWORDS")
\keyword{Monte Carlo summary of aggregated data}% __ONLY ONE__ keyword per line
