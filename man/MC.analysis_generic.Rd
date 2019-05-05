\name{MC.analysis_generic}
\alias{MC.analysis_generic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Analysis of the Monte Carlo simulation (general function)
}
\description{
  General function for running the analysis of the Monte Carlo simulation.
}
\usage{
  MC.analysis_generic(x, delta, qUpper, data.det, sim.det, event.ini, event.end, 
              ntick, summ.data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A \code{list} of 1, which contains the output of the Monte Carlo simulation as 
    a data.frame with n rows as time steps and the first column is time in format POSIXct and m columns named 1, 2, 3... m, where m is the number of Monte Carlo runs results.}
  \item{delta}{
    A \code{numeric} value that specifies the level of aggregation required in minutes.}
  \item{qUpper}{
    A \code{character} string that defines the upper percentile to plot the confidence band of results,
    several options are possible \code{"q999"} the 99.9th percentile, \code{"q995"} the 99.5th percentile,
    \code{"q99"} the 99th percentile, \code{"q95"} the 95th percentile, \code{"q50"} the 50th percentile.
    The lower boundary of the confidence band (showed in gray in the output plots) is the 5th percentile 
    in all cases.}
  \item{data.det}{
    A \code{data.frame} that contains the time series of the 
    main driving force of the system to be simulated deterministically, e.g. precipitation. This 
    data.frame should have 
    only two columns: the first one, Time [y-m-d h:m:s] in POSIXct format; the second one, a numeric value equal to the magnitude
    of the variable.}
  \item{sim.det}{
    A \code{list} of 1 that contains the results of the deterministic simulation, here the output given
    \code{data.det}. The format is the same as \code{data.det}}.
  \item{event.ini}{
    A time-date string in \code{POSIXct} format that defines the initial time for event analysis. }
  \item{event.end}{
    A time-date string in \code{POSIXct} format that defines the final time for event analysis.}
  \item{ntick}{
    A \code{numeric} value to specify the number of ticks in the x-axis for the event time-window plots.}
  \item{summ.data}{
    A \code{list} by default NULL. If provided, the list should contain an output of the MC.analysis function,   
    and the analysis is done again without the calculation of some of the internal variables, therefore the 
    analysis is faster.}
}
%\details{
  %%  ~~ If necessary, more details than the description above ~~
    %}
\value{
  A \code{list} of length 2:
  \item{summ }{A \code{list} that contains the summary statistics of the Monte Carlo simulation per output
  variable. Each output variable is summarised by calculating the mean \code{"Mean"}, standard deviation
  \code{"sd"}, variance \code{"Variance"}, 5th, 25th, 50th, 75th, 95th, 99.5th, 99.9th percentiles 
  \code{"q05", "q25", "q50", "q75", "q95", "q995"}, \code{"q999"}, the max \code{"Max"}, the sum \code{"Sum"}, 
  time \code{"time"}, and the deterministic data \code{"p1"}, all variables as time series.}
  \item{variance }{A \code{data.frame} that contains the summary statistics of the variance of the Monte 
  Carlo simulation per output variable.}
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
  
  \seealso{
    See also \code{\link{setup-class}}, \code{\link{MC.setup-methods}}, \code{\link{MC.sim-methods}}.
  }
\examples{
## Creating meta-model
Model <- function(A, B, variable.1, variable.2){
  lum <- A*variable.1 + B*variable.2
}

## Model input and parameter set-up

time <- data.frame(time = seq.POSIXt(from = as.POSIXct("2019-01-01"), 
                                     to = as.POSIXct("2019-01-02"), by = 60*60*6))
data <- cbind(time, data = 25) 
data

new.setup <- setup(id = "MC_1",
                   nsim = 10,
                   seed = 123,
                   mcCores = 1,
                   ts.input = data, 
                   rng = rng <- list(
                     A = 1.25,
                     B = 0.75,
                     variable.1 = c(pdf = "uni", min = 0, max = 4),
                     variable.2 = c(pdf = "uni", min = 2.2, max = 3.2)
                   )
)

str(new.setup)

## Monte Carlo simulation set-up
set.seed(slot(new.setup, "seed"))

new.mc.setup <- MC.setup(new.setup)
str(new.mc.setup)         

## Monte Carlo simulation
output <- data.frame(time = new.mc.setup$ts.input[,1])
output[,2:(new.mc.setup$nsim + 1)] <- NA

for(i in 1:new.mc.setup$nsim){
  for(j in 1:nrow(new.mc.setup$ts.input)){
    
    ## model parameter definition
    A <- new.mc.setup$par$A
    B <- new.mc.setup$par$B
    
    ## model input definition
    variable.1 <- new.mc.setup$par$variable.1[i,j]
    variable.2 <- new.mc.setup$par$variable.2[i,j]
    
    ## model evaluation
    output[j,i+1] <- Model(A, B, variable.1, variable.2)  
  }
}

output <- list(output1 = output)
output  

## Deterministic simulation
# model parameter definition
A <- new.mc.setup$par$A
B <- new.mc.setup$par$B

# model input definition
variable.1.det <- apply(X = new.mc.setup$par$variable.1, MARGIN = 2, FUN = mean)
variable.2.det <- apply(X = new.mc.setup$par$variable.2, MARGIN = 2, FUN = mean)

output.det      <- Model(A, B, variable.1.det, variable.2.det)  
output.det      <- cbind(time, output.det)
output.det      <- list(out1 = output.det)
str(output.det)

## Monte Carlo analysis
delta     <- 60*6 # minutes
qUpper    <- "q95"
event.ini <- data$time[1]
event.end <- data$time[nrow(data)]
ntick     <- 1

analysis <- MC.analysis_generic(x = output, delta = delta, qUpper = qUpper, data.det = data,
                                sim.det = output.det, event.ini = event.ini, event.end = event.end,
                                ntick = ntick)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Monte Carlo simulation }
\keyword{ MC.analysis_generic }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line