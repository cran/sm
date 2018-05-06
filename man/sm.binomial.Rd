\name{sm.binomial}
\alias{sm.binomial}
\title{
Nonparametric logistic regression
}
\description{
This function estimates the regression curve using the local likelihood
approach for a vector of binomial observations and an associated vector
of covariate values.
}
\usage{
sm.binomial(x, y, N = rep(1, length(y)), h, \dots)
}
\arguments{
  \item{x}{
    vector of the covariate values
  }
  \item{y}{
    vector of the response values; they must be
    nonnegative integers not larger than those of \code{N}.
  }
  \item{h}{
    the smoothing parameter; it must be positive.
  }
  \item{N}{
    a vector containing the binomial denominators.
    If missing, it is assumed to contain all 1's.
  }
  \item{\dots}{   
      other optional parameters are passed to the \code{sm.options}
  function, through a mechanism which limits their effect only to this
  call of the function; those relevant for this function are the following:
    \code{add}, 
    \code{col}, 
    \code{display}, 
    \code{eval.points}, 
    \code{nbins}, 
    \code{ngrid}, 
    \code{pch}, 
    \code{xlab}, 
    \code{ylab};
  see the documentation of  \code{\link{sm.options}} for their description.
}}
\value{
A list containing vectors with the evaluation points, the corresponding
probability estimates, the linear predictors, the upper and lower points
of the variability bands (on the probability scale) and the standard
errors on the linear predictor scale.
}
\section{Side Effects}{
graphical output will be produced, depending on the value of the
\code{display} parameter.
}
\details{
see Sections 3.4 and 5.4 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis:}
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
  \code{\link{sm.binomial.bootstrap}}, \code{\link{sm.poisson}},
  \code{\link{sm.options}}, \code{\link{glm}}, \code{\link{binning}}
}
\examples{\dontrun{
# the next example assumes that all binomial denominators are 1's
sm.binomial(dose, failure, h=0.5)
# in the next example, (some of) the dose levels are replicated 
sm.binomial(dose, failure, n.trials, h=0.5)
}

with(birth, {
   sm.binomial(Lwt[Smoke=="S"], Low[Smoke=="S"], h=20,
           xlab='mother weight[Smoke=="S"]')
   x<- seq(0,1,length=30)
   y<- rbinom(30,10,prob=2*sin(x)/(1+x))
   sm.binomial(x,y,N=rep(10,30), h=0.25)
})
}
\keyword{nonparametric}
\keyword{smooth}
\keyword{models}
% Converted by Sd2Rd version 1.15.


