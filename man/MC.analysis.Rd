\name{MC.analysis}
\alias{MC.analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Analysis of the Monte Carlo simulation
}
\description{
  Function for running the analysis of the Monte Carlo simulation.
}
\usage{
  MC.analysis(x, delta, qUpper, p1.det, sim.det, event.ini, event.end, 
  ntick, summ.data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A \code{list} .}
  \item{delta}{
    A \code{numeric} value that specifies the level of aggregation required in minutes.}
  \item{qUpper}{
    A \code{character} string that defines the upper percentile to plot the confidence band of results,
    several options are possible \code{"q999"} the 99.9th percentile, \code{"q995"} the 99.5th percentile,
    \code{"q99"} the 99th percentile, \code{"q95"} the 95th percentile, \code{"q50"} the 50th percentile.
    The lower boundary of the confidence band (showed in gray in the output plots) is the 5th percentile 
    in all cases.}
  \item{p1.det}{
    A \code{data.frame} that contains the time series of the 
    main driving force of the system to be simulated deterministically, e.g. precipitation. This 
    data.frame should have 
    only two columns: the first one, Time [y-m-d h:m:s]; the second one, a numeric value equal to the magnitude
    of the environmental variable.}
  \item{sim.det}{
    A \code{list} that contains the results of the deterministic simulation, here the output of EmiStatR given
    \code{p1.det}. See the method \code{EmiStatR} from the homonym package for details.}
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
  time \code{"time"}, and the deterministic precipitation \code{"p1"}, all variables as time series.}
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
## the Monte Carlo simulation: MC.sim
library(EmiStatR)

# library(xts)
# data(Esch_Sure2010)
# P <- IsReg(Esch_Sure2010, format = "%Y-%m-%d %H:%M:%S", tz = "CET")
# P1 <- P[[2]]
# P1 <- P1["2010-08",][1:55]
# P1 <- cbind.data.frame(time=index(P1), P1 = coredata(P1))

data(P1)
P1 <- P1[165:(110*2),]
plot(P1[,2], typ="l")

library(stUPscales)

setting_EmiStatR <-  setup(id       = "MC_sim1",
                           nsim     = 3, # # use a larger number to have 
                                           # a proper confidence band of simulatios
                           seed     = 123, 
                           mcCores  = 1, 
                           ts.input = P1,  
                           rng      = rng <- list(
                             qs   = 150,    # [l/PE/d]
                             CODs = c(pdf = "nor", mu = 4.378, sigma = 0.751),    # log[g/PE/d]
                             NH4s = c(pdf = "nor", mu = 1.473, sigma = 0.410),    # log[g/PE/d]
                             qf   = 0.04,     # [l/s/ha]
                             CODf = 0,              # [g/PE/d]
                             NH4f = 0,              # [g/PE/d]
                             CODr = c(pdf = "nor", mu = 3.60, sigma = 1.45),     # 71 log[mg/l]
                             NH4r = 1,              # [mg/l]
                             nameCSO = "E1",        # [-]
                             id      = 1,           # [-]
                             ns      = "FBH Goesdorf",  # [-]
                             nm      = "Goesdorf",  # [-]
                             nc      = "Obersauer", # [-]
                             numc    = 1,          # [-]
                             use     = "R/I",  # [-]
                             Atotal  = 36,              # [ha]
                             Aimp    = c(pdf = "uni", min = 4.5, max = 25),       # [ha]
                             Cimp    = c(pdf = "uni", min = 0.25, max = 0.95),  # [-]
                             Cper    = c(pdf = "uni", min = 0.05, max = 0.60),  # [-]
                             tfS     = 1,               # [time steps]
                             pe      = 650,             # [PE]
                             Qd      = 5,               # [l/s]
                             Dd      = 0.150,           # [m]
                             Cd      = 0.18,            # [-]
                             V       = 190,             #  [m3]
                             lev.ini = 0.10,            # [m]
                             lev2vol = list(lev = c(.06, 1.10, 1.30, 3.30),   # [m]
                                            vol = c(0, 31, 45, 190))          # [m3]
                           ),
                           ar.model  = ar.model <- list(
                             CODs    = 0.5,         
                             NH4s    = 0.5,
                             CODr    = 0.7),
                           var.model = var.model <- list(
                             inp     = c("", ""), # c("CODs", "NH4s"), # c("", ""),
                             w       = c(0.04778205, 0.02079010),
                             A       = matrix(c(9.916452e-01, -8.755558e-05, 
                                                -0.003189094, 0.994553910), nrow=2, ncol=2),
                             C       = matrix(c(0.009126591, 0.002237936, 
                                                0.002237936, 0.001850941), nrow=2, ncol=2)))

MC_setup <- MC.setup(setting_EmiStatR)

sims <- MC.sim(x = MC_setup, EmiStatR.cores = 0)


## Monte Carlo simulation analysis: MC.analysis

# Deterministic simulation
# Definition of structure 1, E1:

E1 <- list(id = 1, ns = "FBH Goesdorf", nm = "Goesdorf", nc = "Obersauer", numc = 1,
           use = "R/I", Atotal = 36, Aimp = 25.2, Cimp = 0.80, Cper = 0.30,
           tfS = 0, pe = 650, Qd = 5, 
           Dd = 0.150, Cd = 0.18, V = 190, lev.ini = 0.10,
           lev2vol = list(lev = c(.06, 1.10, 1.30, 3.30),   
                          vol = c(0, 31, 45, 190))
           )

# Defining deterministic input:
library(EmiStatR)
# data(P1)

input.det <- input(spatial = 0, zero = 1e-5, 
                    folder = system.file("shiny", package = "EmiStatR"),
                    cores = 0,
                    ww = list(qs = 150, CODs = 104, NH4s = 4.7),
                    inf = list(qf= 0.04, CODf = 0, NH4f = 0),
                    rw = list(CODr = 71, NH4r = 1, stat = "Dahl"),
                    P1 = P1, st = list(E1=E1), export = 0)

# Invoking `EmiStatR` with the deterministic input: 
sim.det   <- EmiStatR(input.det)

# further arguments
delta <- 10 # minutes 
qUpper <- "q999" 

event.ini <- as.POSIXct("2016-01-02 03:20:00")
event.end <- as.POSIXct("2016-01-02 12:30:00")

# uncomment to run:
# new_analysis <- MC.analysis(x = sims, delta = delta, qUpper = qUpper,  p1.det = P1, 
#                             sim.det = sim.det, event.ini = event.ini, event.end = event.end, 
#                             ntick = 5, summ.data = NULL)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Monte Carlo simulation }
\keyword{ MC.analysis }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line