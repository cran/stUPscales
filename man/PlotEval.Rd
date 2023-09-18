\name{PlotEval}
\alias{PlotEval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to execute evaluation plot
}
\description{
This function creates an evaluation plot for the Monte Carlo simulation result.
}
\usage{
PlotEval(eval, ts, gof1, namePlot, pos1, pos2, pos3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{eval}{
A \code{data.frame} with	n observations of  seven variables:
 1) \code{time}: A \code{POSIXct} object with format "\%Y-\%m-\%d \%H:\%M:\%S" defining the time vector;
 2) \code{column 2}: a numeric vector containing the values of the observed variable, which is the first variable of the \code{Level2Volume} relationship;
 3) \code{column 3}: a numeric vector containing the values for the second variable of the \code{Level2Volume} relationship;
 4) \code{column 4}: a numeric vector containing the corresponding simulated values for the second variable of the \code{Level2Volume} relationship;
 5) \code{column 5}: a numeric vector containing the difference between the vectors \code{volT_sim} and \code{volT_obs}.
 6) \code{Rainfall}: a numeric vector named "Rainfall" containing the values of the driving force variable used in the simulations, e.g. rainfall.
 7) \code{column 7}: (Optional) a numeric vector containing the values of the driving force variable used in the simulations in other measurement units, e.g. rainfall in intensity units if rainfall is the driving force of the simulations.
 }
  \item{ts}{
An \code{xts} object representing the \code{eval} \code{data.frame} indexed by the \code{time} vector of the \code{eval} argument:
containing six data variables as it is defined by the \code{eval} argument: 1) \code{column 2}; 2) \code{column 3}; 3) \code{column 4}; 4) \code{column 5}; 5) \code{Rainfall}; 6) \code{column 7}.
}
  \item{gof1}{
A matrix with the output of \code{GoF} function.
}
  \item{namePlot}{
A character string defining the name of the plot to be created.
}
  \item{pos1}{
Location to place the legend on the inside of the first sub-plot frame. Can be one of "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". 
}
  \item{pos2}{
Location to place the legend on the inside of the second sub-plot frame. Can be one of "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". 
}
  \item{pos3}{
Location to place the legend on the inside of the third sub-plot frame. Can be one of "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The function creates a plot in the current working directory with the goodness-of-fit between simulations and observations. 
The plot is provided in pdf format. 
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
time <- seq(from = as.POSIXct("2017-11-09"), by = 60*60*24, length.out = 230) # the time vector
data <- cbind.data.frame(time, NA) # a NA vector
data[,3] <- rnorm(230, .25, .1) # random normal distributed data, obs
data[,4] <- data[,3]*1.2  # positive correlated data, sim
data[,5] <- data[,4] - data[,3] # difference sim and obs
data[,6] <- 0 # driving force
data[,7] <- NA # a NA vector

colnames(data) <- c("time", "var1", "obs", "sim", "difference", "Rainfall", "Rainfall2")
head(data)

ts <- IsReg.ts(data, "\%Y-\%m-\%d", "ECT")
ts <- ts[[2]]

gof.new <- GoF(data, 4, 3, "")
gof.new

## not run
## creating the plot (uncomment to run)
#PlotEval(data, ts, gof.new, "ExamplePlot", "topright", "topright", "topright")
#}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{PlotEval}% use one of  RShowDoc("KEYWORDS")
\keyword{Plot}