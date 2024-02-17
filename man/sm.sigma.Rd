\name{sm.sigma}
\alias{sm.sigma}
\title{Estimation of the error standard deviation in nonparametric regression.}
\description{This function estimates the error standard deviation in nonparametric
       regression with one or two covariates.}
\usage{sm.sigma(x, y, rawdata = NA, weights = rep(1, length(y)), 
               diff.ord = 2, ci = FALSE, model = "none", h = NA, \dots)
       }
\arguments{
  \item{x}{a vector or two-column matrix of covariate values.}
  \item{y}{a vector of responses.}
  \item{rawdata}{a list containing the output from a binning operation.
           This argument is used by \code{sm.regression} and it need not be
           set for direct calls of the function.}
  \item{weights}{a list of frequencies associated with binned data.
           This argument is used by \code{sm.regression} and it need not be
           set for direct calls of the function.}
  \item{diff.ord}{an integer value which determines first (1) or second (2)
           differencing in the estimation of sigma.}
  \item{ci}{a logical value which controls whether a confidence interval is
           produced.}
  \item{model}{a character variable.  If this is set to \code{"constant"}
           then a test of constant variance over the covariates is performed
           (only in the case of two covariates)}
  \item{h}{a vector of length two defining a smoothing parameter to be used
           in the test of constant variance.}
  \item{\dots}{other optional parameters are passed to the \code{sm.options}
           function, through a mechanism which limits their effect only to this
           call of the function; the only one relevant for this function is
           \code{nbins}.}
  }
\value{a list containing the estimate and, in the two covariate case, a 
           matrix which can be used by the function \code{sm.sigma2.compare}, 
           pseudo-residuals and, if appropriate, a confidence interval and 
           a p-value for the test of constant variance.}
\section{Side Effects}{none.}
\details{see the reference below.}
\references{Bock, M., Bowman, A.W. & Ismail, B. (2007).
           Estimation and inference for error variance in bivariate 
           nonparametric regression.
           \emph{Statistics & Computing}, to appear.}
\seealso{\code{\link{sm.sigma2.compare}}}
\examples{
\dontrun{
with(airquality, {
   x     <- cbind(Wind, Temp)
   y     <- Ozone^(1/3)
   group <- (Solar.R < 200)
   sig1 <- sm.sigma(x[ group, ], y[ group], ci = TRUE)
   sig2 <- sm.sigma(x[!group, ], y[!group], ci = TRUE)
   print(c(sig1$estimate, sig1$ci))
   print(c(sig2$estimate, sig2$ci))
   print(sm.sigma(x[ group, ], y[ group], model = "constant", h = c(3, 5))$p)
   print(sm.sigma(x[!group, ], y[!group], model = "constant", h = c(3, 5))$p)
   print(sm.sigma2.compare(x[group, ], y[group], x[!group, ], y[!group]))
})
}}
\keyword{nonparametric}
\keyword{smooth}
