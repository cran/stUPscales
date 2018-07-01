\name{GoF}
\alias{GoF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper function for the \code{gof} function from \code{hydroGOF} package
}
\description{
A wrapper function for the \code{gof} function from \code{hydroGOF} package
}
\usage{
GoF(eval, col_sim, col_obs, name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eval}{
A \code{matrix} or \code{data.frame} with	n observations of at least two variables: simulations and observations.
 }
  \item{col_sim}{
A numeric value defining the column in \code{eval} data.frame that contains the simulated vector time series.
}
  \item{col_obs}{
A numeric value defining the column in \code{eval} data.frame that contains the observed vector time series.
}
  \item{name}{
A character string that defines the name of the files (.csv and .RData) created with the results. If missing then no files are created.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A vector with 20 elements for each one of the following measures of godness-of-fit: 1) ME, mean error; 2) MAE, mean absolute error; 
3) MSE, mean squared error; 4) RMSE, root mean square error; 5) NRMSE \%, normalized root square error; 6) PBIAS \%, percent bias; 
7) RSR, Ratio of RMSE to the standard deviation of the observations; 8) rSD, Ratio of Standard Deviations; 
9) NSE, Nash-Sutcliffe Efficiency; 10) mNSE, modified Nash-Sutcliffe efficiency; 11) rNSE, relative Nash-Sutcliffe efficiency; 
12) d, Index of Agreement; 13) md, Modified index of agreement; 14) rd, Relative Index of Agreement; 15) cp,Coefficient of persistence; 
16) r, Pearson product-moment correlation coefficient; 17) R2, Coefficient of Determination; 
18) bR2, Coefficient of determination (r2) multiplied by the slope of the regression line between sim and obs; 
19) KGE,Kling-Gupta Efficiency; 20) VE, Volumetric Efficiency. 
}
\references{
Mauricio Zambrano-Bigiarini, 2014. hydroGOF: Goodness-of-fit functions for comparison of simulated and observed hydrological time series.
R package version 0.3-8. \cr
https://CRAN.R-project.org/package=hydroGOF.
}
\author{
J.A. Torres-Matallana}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
library(stUPscales)

data_new <- rnorm(230, .25, .1)  
data_new <- cbind(data_new, data_new*1.2)
colnames(data_new) <- c("sim", "obs")
head(data_new)

gof.new <- GoF(data_new, 1, 2)
gof.new

# writing files
gof.new <- GoF(data_new, 1, 2, "GoF_results")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{GoF}% use one of  RShowDoc("KEYWORDS")
\keyword{goodness-of-fit}% __ONLY ONE__ keyword per line
\keyword{hydroGOF}% __ONLY ONE__ keyword per line