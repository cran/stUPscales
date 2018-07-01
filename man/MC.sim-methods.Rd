\name{MC.sim-methods}
\docType{methods}
\alias{MC.sim-methods}
\alias{MC.sim,list,numeric-method}
\alias{MC.sim}
\title{ ~~ Methods for Function \code{MC.sim}  ~~}
\description{
Method to be invoked for running the Monte Carlo simulation. The simulator used is the method \code{EmiStatR}
from the homonym package. This method should be rewritted for working with another simulator.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "list", EmiStatR.cores = "numeric")}}{
%%  ~~describe this method here~~
}
}}

\usage{
MC.sim(x, EmiStatR.cores)
}

\arguments{
  \item{x}{
   an object of class \code{list} as is defined by method \code{MC.setup}.}
  \item{EmiStatR.cores}{
   a \code{numeric} value for specifying the number of cores (CPUs) to be used in the EmiStatR method. 
   Use zero for not use parallel computation. See class \code{input} of package \code{EmiStatR} for details. }
}

\value{
  A \code{list} of length 2:
  \item{mc }{A \code{list} that contains the \code{MC_setup}, \code{timing} and \code{lap} objects.}
  \item{sim1 }{A \code{list} that contains the Monte Carlo matrices of the simulator output.}
}

\examples{
## the Monte Carlo simulation: MC.sim

library(EmiStatR)
data(P1)
P1 <- P1[165:(110*2),]
plot(P1[,2], typ="l")

library(stUPscales)

setting_EmiStatR <-  setup(id       = "MC_sim1",
                           nsim     = 3, # use a larger number to have 
                                         # a proper confidence band of simulations
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
str(sims)
}

\keyword{methods}
\keyword{Monte Carlo simulation}
\keyword{MC.sim}
