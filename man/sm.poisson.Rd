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
sm.poisson(x, y, h, ngrid=25, eval.points, add=F, 
           display="estimate", xlab, ylab, pch=1, col=2, \dots)
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
\item{ngrid}{
the number of points where the regression curve must be estimated
(only used if \code{eval.points} is not given).
}
\item{eval.points}{
the vector of points on the \code{x} axis where the regression must be
estimated. If the parameter \code{eval.points} is not given, this vector
is chosen to be formed by \code{ngrid} equally spaced points between
\code{min(x)} and \code{max(x)}.
}
\item{add}{
if graphical output is produced, this parameter controls whether a 
new plot is created, or graphical output is added to the existing one.
}
\item{display}{
controls the type of graphical output; possible values are
\code{"estimate"} (default), \code{"se"}, \code{"none"}.
}
\item{xlab}{
label of the x-axis
}
\item{ylab}{
label of the y-axis
}
\item{pch}{
plotting character of the raw observed frequency.
}
\item{col}{
colour used for plotting curves and points
}
\item{...}{
additional graphical parameters
}}
\value{
A list containing vectors with the evaluation points, the corresponding
probability estimates, the linear predictors, the upper and lower points
of the variability bands  and the standard errors on the linear predictor 
scale.
}
\section{Side Effects}{
graphical output will be produced, depending on the value of the
display parameter, unless this is set to \code{"none"}.
}
\details{
see Sections 3.4 and 5.4 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.logit}}, \code{\link{sm.logit.bootstrap}}, \code{\link{sm.poisson.bootstrap}}
}
\examples{
sm.poisson(exposure.time, N.events, 0.5, display="se")
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 0.3-2.
