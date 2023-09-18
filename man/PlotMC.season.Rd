\name{PlotMC.season}
\alias{PlotMC.season}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A plot function for time series seasons
}
\description{
This is an internal function invoked by \code{MC.analysis} function to generate a season plot of the time series under analysis. 
A season means a time series with length greater to one month e.g. montly, yearly, decadal time series.
}
\usage{
PlotMC.season(summ1, namePlot, ylab, qUpper)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summ1}{
A data.frame with n observations of m variables as is provided by the output of function \code{MC.summary.agg}
for the variable to be plotted, which the summary was computed.}
  \item{namePlot}{
A character string defining the name of the plot. The file created with the plot has this name.
}
  \item{ylab}{
A character string to define the label of the axes y for the variable to plot.
}
  \item{qUpper}{
A character string that defines the upper percentile to plot the confidence band
of results, several options are possible "q999" the 99.9th percentile, "q995" the
99.5th percentile, "q99" the 99th percentile, "q95" the 95th percentile, "q50"
the 50th percentile. The lower boundary of the confidence band (showed in gray
in the output plots) is the 5th percentile in all cases.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The function creates the plot in the current working directory. The format of the plot is pdf.
}

%\references{
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

data("P1")
P1 <- P1[1:550,]
new_data <- matrix(data = NA, nrow = nrow(P1), ncol = 55)
for(i in 1:55){
  new_data[,i] <- matrix(data = rnorm(nrow(P1), 22, 11), 
                         nrow = nrow(P1), ncol = 1)
}
new_data <- t(new_data)

new_summary <- MC.summary(p1 = P1, data = new_data)
head(new_summary)
dim(new_summary)

new_summary$month <- strftime(new_summary[,"time"], format = "\%Y-\%m")

## not run
## creating the plot (uncomment to run)
# PlotMC.season(summ1 = new_summary, namePlot = "ExamplePlot", 
#              ylab = "Variable 1 [units]", qUpper = "q95")
#}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{PlotMC.season}% use one of  RShowDoc("KEYWORDS")
\keyword{Monte Carlo simulation}% __ONLY ONE__ keyword per line
\keyword{Plot}% __ONLY ONE__ keyword per line
