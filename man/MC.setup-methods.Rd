\name{MC.setup-methods}
\docType{methods}
\alias{MC.setup-methods}
\alias{MC.setup,setup-method}
\alias{MC.setup}
\title{Methods for Function \code{MC.setup}}
\description{
Given an object of class \code{setup}, 
the method can be invoked for setting-up the Monte Carlo simulation. 
The variables are sampled accordingly to their parameters specified in the slot \code{rng} of the 
\code{setup} object. If \code{ar.model} is defined in slot \code{ar.model}, then the specified
variables are sampled from the \code{pdf} \code{nor} as an autorregresive (AR) model via the
function \code{arima.sim} from base package \code{stats}. If \code{var.model} is defined in 
slot \code{var.model}, then the specified variables are sampled from the \code{pdf} \code{nor} as 
an vector autorregresive (VAR) model via the function \code{mAr.sim} from package \code{mAr} 
(see Barbosa, 2015, and Luetkepohl, 2005, for details). See \code{setup-class} for further details 
to define the AR and VAR models.
}

\section{Methods}{
\describe{

\item{\code{signature(x = "setup")}}{


}
}}

\usage{
MC.setup(x)
         
}

\arguments{
  \item{x}{
   an object of class \code{setup}.
}}


\author{
J.A Torres-Matallana
}

\references{
S. M. Barbosa, Package "mAr": Multivariate AutoRegressive analysis, 1.1-2, The Comprehensive R Archive Network, CRAN, 2015.

H. Luetkepohl, New Introduction to Multiple Time Series Analysis, Springer, 2005.
}

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

MC_setup <- MC.setup(ini)
str(MC_setup)

## definition of AR models for variables var2 and var3 with AR coefficients 0.995 and 0.460

library(EmiStatR)
data(P1)

ini_ar <- setup(id = "MC_sim1_ar", nsim = 500, seed = 123, mcCores = 1, ts.input = P1,
                rng = list(var1 = 150, var2 = c(pdf = "nor", mu = 150, sigma = 5),
                           var3 = c(pdf = "nor", mu = 90, sigma = 2.25)),
                ar.model = ar.model <- list(var2 = 0.995, var3 = 0.460)
)

MC_setup_ar <- MC.setup(ini_ar)
str(MC_setup_ar)

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


MC_setup_var <- MC.setup(ini_var)
str(MC_setup_var)
}


\keyword{methods}
\keyword{Monte Carlo simulation}
\keyword{MC.setup}
\keyword{setup}