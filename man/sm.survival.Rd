\name{sm.survival}
\alias{sm.survival}
\title{
Nonparametric regression with survival data.
}
\description{
This function creates a smooth, nonparametric estimate of the quantile
of the distribution of survival data as a function of a single covariate.
A weighted Kaplan-Meier survivor function is obtained by smoothing across 
the covariate scale.  A small amount of smoothing is then also applied across 
the survival time scale in order to achieve a smooth estimate of the quantile.
}
\usage{
sm.survival(x, y, status, h , hv = 0.05, p = 0.5, status.code = 1, \dots)
}
\arguments{
\item{x}{
a vector of covariate values.
}
\item{y}{
a vector of survival times.
}
\item{status}{
an indicator of a complete survival time or a censored value.  The value of
\code{status.code} defines a complete survival time.
}
\item{h}{
the smoothing parameter applied to the covariate scale.  A normal kernel
function is used and \code{h} is its standard deviation.
}
\item{hv}{
a smoothing parameter applied to the weighted Kaplan-Meier functions derived
from the smoothing procedure in the covariate scale.  This ensures that
a smooth estimate is obtained.
}
\item{p}{
the quantile to be estimated at each covariate value.
}
\item{status.code}{
the value of \code{status} which defines a complete survival time.
}
\item{...}{
other optional parameters are passed to the \code{sm.options} function, through
a mechanism which limits their effect only to this call of the function; 
those relevant for this function are the following:
}
\item{eval.points}{
the points at which the estimate will be evaluated.
}
\item{ngrid}{
the number of points in a regular grid over the covariate scale at which
the estimate will be evaluated, if \code{eval.points} is set to \code{NA}.
Default: \code{ngrid=50}.
}
\item{display}{
The setting \code{"none"} will prevent any graphical output from being produced.  
The default setting \code{"lines"} (or indeed any other value) will produce a
plot of the data and estimate.
}
\item{xlab}{
the label attached to the x-axis.
}
\item{ylab}{
the label attached to the y-axis.
}
\item{lty}{
the line type used to plot the estimate.
}
\item{add}{
a logical value which controls whether the estimate is added to
the current plot.  Its default value is F, which creates a new plot.
}}
\value{
a list containing the values of the estimate at the evaluation points
and the values of the smoothing parameters for the covariate and survival
time scales.
}
\section{Side Effects}{
none.
}
\details{
see Section 3.5 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.regression}}, \code{\link{sm.options}}
}
\examples{
x <- runif(50, 0, 10)
y <- rexp(50, 2)
z <- rexp(50, 1)
status <- rep(1, 50)
status[z<y] <- 0
y <- pmin(z, y)
sm.survival(x, y, status, h=2)
}
\keyword{nonparametric}
\keyword{smooth}
\keyword{survival}
% Converted by Sd2Rd version 1.15.
