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
other optional parameters are passed to the \code{sm.options} function, through
a mechanism which limits their effect only to this call of the function; 
those relevant for this function are the following:
\itemize{
\item{add}{
if graphical output is produced, this parameter controls whether a 
new plot is created, or graphical output is added to the existing one.
}
\item{col}{
colour used for plotting curves and points. Default: \code{col=1}.
}
\item{display}{
controls the type of graphical output; possible values are
\code{"estimate"} (default), \code{"se"}, \code{"none"}.
}
\item{eval.points}{
the vector of points on the \code{x} axis where the regression must be
estimated. If the parameter \code{eval.points} is not given, this vector
is chosen to be formed by \code{ngrid} equally spaced points between
\code{min(x)} and \code{max(x)}.
}
\item{nbins}{
The number of bins used when binning operation is performed.
If \code{nbins=0}, binning is not performed; if \code{nbins=NA} (default),
binning is switched on when the number of design points exceeds 100.
}
\item{ngrid}{
the number of points where the regression curve must be estimated
(only used if \code{eval.points} is not given). Default: \code{ngrid=25}.
}
\item{pch }{
plotting character of the raw observed frequency. Default: \code{pch=1}.
}
\item{xlab}{
label of the x-axis. Default is the name of \code{x} object.
}
\item{ylab}{
label of the y-axis .Default is the name of \code{y} object.
}}}}
\value{
A list containing vectors with the evaluation points, the corresponding
probability estimates, the linear predictors, the upper and lower points
of the variability bands  and the standard errors on the linear predictor 
scale.
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
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.binomial}}, \code{\link{sm.binomial.bootstrap}}, \code{\link{binning}}, \code{\link{glm}}
}
\examples{
\dontrun{sm.poisson(exposure.time, N.events, 0.5, display="se")}
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
