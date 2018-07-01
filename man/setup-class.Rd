\name{setup-class}
\Rdversion{1.1}
\docType{class}
\alias{setup-class}
\alias{setup}
\alias{stUPscales,setup-method}

\title{Class \code{"setup"}}
\description{
Class to create objects of signature setup. setup object should be passed to the method
MC.setup.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{setup()}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{id}:}{Object of class \code{"character"} to identify the Monte Carlo simulation.}
    \item{\code{nsim}:}{Object of class \code{"numeric"} to specify the number of Monte Carlo runs.}
    \item{\code{seed}:}{Object of class \code{"numeric"} to specify the seed of the random numbers 
    generator.}
    \item{\code{mcCores}:}{Object of class \code{"numeric"} to specify the number of cores (CPUs) to be
    used in the Monte Carlo simulation.}
    \item{\code{ts.input}:}{Object of class \code{"data.frame"} that contains the time series of the 
    main driving force of the system to be simulated, e.g. precipitation. This data.frame should have 
    at least two columns: the first one, Time [y-m-d h:m:s]; the second one, a numeric value equal to the 
    magnitude of the environmental variable. This data.frame
    can also contain more that one column to allow several time series in several columns. If the data.frame
    has more than two columns, then the number of columns should be at least equal to nsim. If the
    number of columns is greater than nsim, the columns in excess are not recycled because the simulation
    will last nsim iterations.
    }
    \item{\code{rng}:}{Object of class \code{"list"} that contains the names and values of the variables
    to be used in the Monte Carlo simulation. Five modes are available: 1) constant value, i.e. this variable
    will have a constant value along the Monte Carlo simulation; 2) a variable sampled from a uniform (uni)
    probability distribution function (pdf) with parameters for the lower boundary min and upper boundary max; 
    3) a variable sampled from a normal (nor) pdf with parameteres mean mu and standard deviation sigma;
    4) a variable sampled from an autorregresive  (AR) model and normal (nor) pdf with parameteres mean mu and
    standard deviation sigma, the coefficients of the AR model should be defined in the slot ar.model;
    5) a variable sampled from an vector autorregresive (VAR) model and normal (nor) pdf with parameteres mean
    mu and standard deviation sigma, this mode is enabled by defining the vector of intercept terms w, 
    the matrix of AR coefficients A, and the noise covariance matrix C in the slot var.model. 
    See examples for the definition of this slot.}
    \item{\code{ar.model}:}{Object of class \code{"list"} containing the coefficients of the AR model as vectors 
    which name is the variable to be modeled and length the order of the model as is required for function
    arima.sim from the base package stats. 
    The named variables here should correspond to a pdf nor in the slot rng. See examples for the 
    definition of this slot.}
    \item{\code{var.model}:}{Object of class \code{"list"} containing the the vector of intercept terms w, 
    the matrix of AR coefficients A, and the noise covariance matrix C of the VAR model which name is the 
    variable to be modeled and length the order of the model as is required for function
    mAr.sim from the package mAr. The named variables in this slot should correspond to 
    a pdf nor in the slot rng. The current implementation considers the bi-variate case. See examples for 
    the definition of this slot. For mathemamtical details see Luetkepohl (2005).}
  }
}
\section{Methods}{
  \describe{
    \item{MC.setup}{\code{signature(x = "setup")}: execute {MC.setup} function }
    }
}

\references{
S. M. Barbosa, Package "mAr": Multivariate AutoRegressive analysis, 1.1-2, The Comprehensive R Archive Network, CRAN, 2015.

H.Luetkepohl, New Introduction to Multiple Time Series Analysis, Springer, 2005.
}
\author{
J.A Torres-Matallana
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}

\examples{
# loading a precipitation time series as input for the setup class

library(EmiStatR)
data(P1)

# A setup with three variables to be considered in the Monte Carlo simulation: 
# var1, a constant value variable; var2, a variable sampled from a uniform (uni)
# probability distribution function (pdf) with parameters min and max; 
# var3, a variable sampled from a normal (nor) pdf with parameteres mu and sigma

ini <- setup(id = "MC_sim1", nsim = 500, seed = 123, mcCores = 1, ts.input = P1,
             rng = list(var1 = 150, var2 = c(pdf = "uni", min = 50, max = 110),
                        var3 = c(pdf = "nor", mu = 90, sigma = 2.25))
)

str(ini)

## definition of AR models for variables var2 and var3 with AR coefficients 0.995 and 0.460

library(EmiStatR)
data(P1)

ini_ar <- setup(id = "MC_sim1_ar", nsim = 500, seed = 123, mcCores = 1, ts.input = P1,
                rng = list(var1 = 150, var2 = c(pdf = "nor", mu = 150, sigma = 5),
                           var3 = c(pdf = "nor", mu = 90, sigma = 2.25)),
                ar.model = ar.model <- list(var2 = 0.995, var3 = 0.460)
)

str(ini_ar)

## definition of a bi-variate VAR model for variables var2 and var3

ini_var <- setup(id = "MC_sim1_ar", nsim = 500, seed = 123, mcCores = 1, ts.input = P1,
                 rng = rng <- list(var1 = 150, 
                                   var2 = c(pdf = "nor", mu = 150, sigma = 5),
                                   var3 = c(pdf = "nor", mu = 90, sigma = 2.25)),
                 var.model = var.model <- list( inp = c("var2", "var3"), 
                                                w = c(0.048, 0.021),
                                                A = matrix(c(0.992, -8.8e-05, -31e-4, 0.995), 
                                                nrow=2, ncol=2),
                                                C = matrix(c(0.0091, 0.0022, 0.0022, 0.0019), 
                                                nrow=2, ncol=2))
)

str(ini_var)
}

\keyword{classes}

%\usage{
%setup(id = "MC_sim_1", nsim = 1, seed = 0.7010607, mcCores = 1, ts.input = NULL,
%      rng = NULL, ar.model  = list(NULL), var.model = list(NULL))
%}

%\arguments{
%  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
%}
%}
