\name{sm.ancova}
\alias{sm.ancova}
\title{
Nonparametric analysis of covariance
}
\description{
This function allows a set of nonparametric regression curves to be
compared, both graphically and formally in a hypothesis test.  A reference
model, used to define the null hypothesis, may be either equality or 
parallelism.  Regression surfaces can also be compared in a test but a
graphical display is not produced.  
}
\usage{
sm.ancova(x, y, group, h, model = "none", h.alpha = NA, weights=NA,
                 covar = diag(1/weights), ...)

}
\arguments{
\item{x}{
a vector or two-column matrix of covariate values.
}
\item{y}{
a vector of response values.
}
\item{group}{
a vector of group indicators.
}
\item{h}{
the smoothing parameter to be used in the construction of each of the
regression curves.  If this is missing the method of smoothing parameter
selection specified by \code{sm.options} will be applied.
}
\item{model}{
a character variable which defines the reference model.  The values
\code{"none"}, \code{"equal"} and \code{"parallel"} are possible.
}
\item{h.alpha}{
the value of the smoothing parameter used when estimating the vertical
separations of the curves under the parallelism model.  If this is missing,
it is set to 2 * r / n, where r denotes the range of the data and n the
sample size.
}
\item{weights}{
  case weights; see the documentation of \code{\link{sm.regression}}
  for a full description.
}
\item{covar}{
  the (estimated) covariance matrix of y.  The default value assumes
  the data to be independent.  Where appropriate, the covariance structure
  of \code{y} can be estimated by the user, externally to \code{sm.ancova},
  and passed through this argument.  This is used in the hypothesis tests
  but not in the construction of the reference band for comparing two groups
  (and so graphics are disabled in this case).
}
\item{\dots}{
  other optional parameters are passed to the \code{sm.options}
  function, through a mechanism which limits their effect only to this
  call of the function. Those relevant for this function are the following:
  \code{display},
  \code{ngrid},
  \code{eval.points},
  \code{xlab},
  \code{ylab};
see the documentation of  \code{\link{sm.options}} for their description.
}}
\value{
a list containing an estimate of the error standard deviation and, where
appropriate, a p-value and reference model.  If the parallelism model has
been selected then a vector of estimates of the vertical separations of the 
underlying regression curves is also returned.  If a reference band has been
requested, the upper and lower boundaries and their common evaluation points
are also returned.
}
\section{Side Effects}{
a plot on the current graphical device is produced, unless \code{display="none"}
}
\details{
see Sections 6.4 and 6.5 of the book by Bowman \& Azzalini, and
the papers by Young \& Bowman listed below.
This function is a developed version of code originally written by Stuart Young.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.

Young, S.G. and Bowman, A.W. (1995).  
Nonparametric analysis of covariance.  
\emph{Biometrics}
\bold{51}, 920--931.

Bowman, A.W. and Young, S.G. (1996).  
Graphical comparison of nonparametric curves.  
\emph{Applied Statistics}
\bold{45}, 83--98.
}
\seealso{
\code{\link{sm.regression}}, \code{\link{sm.density.compare}}, \code{\link{sm.options}}
}
\examples{
x <- runif(50, 0, 1)
y <- 4*sin(6*x) + rnorm(50)
g <- rbinom(50, 1, 0.5)
sm.ancova(x, y, g, h = 0.15, model = "equal")
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
