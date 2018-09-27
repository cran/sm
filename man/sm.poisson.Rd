\name{sm.poisson}
\alias{sm.poisson}
\title{
Nonparametric Poisson regression
}
\description{
This function estimates the regression curve using the local likelihood
approach for a vector of Poisson observations and an associated vector
of covariate values.
}
\usage{
sm.poisson(x, y, h, \dots)
}
\arguments{
\item{x}{
vector of the covariate values
}
\item{y}{
vector of the response values; they must be nonnegative integers.
}
\item{h}{
the smoothing parameter; it must be positive.
}
\item{\dots}{
 other optional parameters are passed to the \code{sm.options} function,
through a mechanism which limits their effect only to this call of the
function. Those specifically relevant for this function are the following:
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
of the variability bands  and the standard errors on the linear predictor 
scale.
}
\section{Side Effects}{
graphical output will be produced, depending on the value of the
\code{display} option.
}
\details{
see Sections 3.4 and 5.4 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997).
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
  \code{\link{sm.binomial}}, \code{\link{sm.binomial.bootstrap}},
  \code{\link{binning}}, \code{\link{glm}}
}
\examples{
with(muscle, {
   TypeI <- TypeI.R+ TypeI.P+TypeI.B
   sm.poisson(x=log(TypeI), y=TypeII, h=0.25,display="se")
   sm.poisson(x=log(TypeI), y=TypeII, h=0.75, col=2, add=TRUE)
})
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
