\name{PlotMC.event}
\alias{PlotMC.event}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A plot function for time series events}
\description{
This is an internal function invoked by \code{MC.analysis} function to generate an event plot of the time series under analysis.
A event means a time series with length lower to one month i.e. sub-montly time series.
}
\usage{
PlotMC.event(summ, summ1, obs, det.var, det.var1, namePlot, ylab, ylab1, ntick, qUpper)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summ}{
A data.frame with n observations of m variables as is provided by the output of function \code{MC.summary.agg}
for the first variable to be plotted.
}
  \item{summ1}{
A data.frame with n observations of m variables as is provided by the output of function \code{MC.summary.agg}
for the second variable to be plotted.
}
  \item{obs}{
A numeric value equal to 0. used for internal use.
}
  \item{det.var}{
A character string defining the name of the first variable from \code{summ} object to be plotted.
}
  \item{det.var1}{
A character string defining the name of the second variable from \code{summ} object to be plotted.
}
  \item{namePlot}{
A character string defining the name of the plot. The file created with the plot has this name.
}
  \item{ylab}{
A character string to define the label of the axes y for the first variable sub-plot.
}
  \item{ylab1}{
A character string to define the label of the axes y for the second variable sub-plot.
}
  \item{ntick}{
A numeric value integer which defines the number of tick marks in the axis x of the sub-plots.
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

# definition of the first summary.agg object
data("P1")
P1 <- P1[1:1100,]

new_data <- matrix(data = NA, nrow = nrow(P1), ncol = 55)
for(i in 1:55){
  new_data[,i] <- matrix(data = rnorm(nrow(P1), 45, 15), 
                         nrow = nrow(P1), ncol = 1)
}
new_data <- t(new_data)

new_summary <- MC.summary(p1 = P1, data = new_data)

# deterministic simulation
det <- rnorm(nrow(P1), 45, 15)
det <- cbind(det, rnorm(nrow(P1), 55, 23))
colnames(det) <- c("det1", "det2")

# level of aggregation
delta <- 60*2 # 2 hours

new_summary_agg <- MC.summary.agg(summ = new_summary, det, delta,
                                  func.agg = mean, func.agg.p = sum)

# definition of the second summary.agg object
new_data1 <- matrix(data = NA, nrow = nrow(P1), ncol = 55)
for(i in 1:55){
  new_data1[,i] <- matrix(data = rnorm(nrow(P1), 55, 23), 
                         nrow = nrow(P1), ncol = 1)
}
new_data1 <- t(new_data1)

new_summary1 <- MC.summary(p1 = P1, data = new_data1)
new_summary_agg1 <- MC.summary.agg(summ = new_summary1, det, delta, 
                                   func.agg = mean, func.agg.p = sum)

## not run
## creating the plot for the event (uncomment to run)
#PlotMC.event(summ = new_summary_agg, summ1 = new_summary_agg1, obs = 0, 
#             det.var = "det1", det.var1 = "det2", namePlot = "ExamplePlot", 
#             ylab = "Variable 1 [units]", ylab1 = "Variable 2 [units]", 
#             ntick=10, qUpper= "q95")
#}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{PlotMC.event}% use one of  RShowDoc("KEYWORDS")
\keyword{Monte Carlo simulation}% __ONLY ONE__ keyword per line
\keyword{Plot}% __ONLY ONE__ keyword per line