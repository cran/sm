\name{sm.sigma2.compare}

\alias{sm.sigma2.compare}

\title{Comparison across two groups of the error standard deviation in 
       nonparametric regression with two covariates.}
       
\description{This function compares across two groups, in a hypothesis
       test, the error standard deviation in nonparametric regression 
       with two covariates.}
       
\usage{sm.sigma2.compare(x1, y1, x2, y2)}

\arguments{
  \item{x1}{a two-column matrix of covariate values for group 1.}
  \item{y1}{a vector of responses for group 1.}
  \item{x2}{a two-column matrix of covariate values for group 2.}
  \item{y2}{a vector of responses for group 2.}
  }
  
\value{a p-value for the test of equality of standard deviations.}

\section{Side Effects}{none.}

\details{see the reference below.}

\references{Bock, M., Bowman, A.W.\ \& Ismail, B. (2007).
           Estimation and inference for error variance in bivariate 
           nonparametric regression.
           \emph{Statistics \& Computing}, to appear.}
           
\seealso{\code{\link{sm.sigma}}}

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
