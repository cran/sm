\name{sm.regression.autocor}
\alias{sm.regression.autocor}
\title{
Nonparametric regression with autocorrelated errors
}
\description{
This function estimates nonparametrically the regression function
of \code{y} on \code{x} when the error terms are serially correlated.
}
\usage{
sm.regression.autocor(x=1:n, y, h.first, minh, maxh, method="direct",
                      ngrid=15, display="plot", \dots)
}
\arguments{
\item{y}{
vector of the reponse values
}
\item{h.first}{
the smoothing parameter used for the initial smoothing stage.
}
\item{x}{
vector of the covariate values; if unset, it is assumed to be \code{1:length(y)}.
}
\item{minh}{
the minimum value of the interval where the optimal smoothing parameter
is searched for (default is 0.5).
}
\item{maxh}{
the maximum value of the interval where the optimal smoothing parameter
is searched for (default is 10).
}
\item{method}{
character value which specifies the optimality criterium adopted;
possible values are \code{"no.cor"}, \code{"direct"} (default), and \code{"indirect"}.
}
\item{ngrid}{
the number of points to be considered in \code{(minh,maxh)}.
}
\item{display}{
if this is equal to \code{"plot"}, graphical output is produced on the
current graphical device.
}
\item{...}{
additional graphical parameters.
}}
\value{
a list as returned from sm.regression called with the new value of
smoothing parameter, with an additional term \code{$aux} added which contains
the initial value \code{h.first}, the estimated curve using \code{h.first}, 
the autocorrelation function of the residuals from the initial fit, 
and the residuals.
}
\section{Side Effects}{
a new suggested value for h is printed, and, if the parameter display
is equal to \code{"plot"}, graphical output is produced on the current graphical
device.
}
\details{
see Section 7.5 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.regression}}, \code{\link{sm.autoregression}}
}
\examples{
\dontrun{
	type("Suggested value of h: ", h1)
	sm1 <- sm.regression.eval.1d(x, y, h = h1, hmult = 1, model = "none")
	if(missing(x))
		x.name <- "time"
	else x.name <- deparse(substitute(x))
	if(display == "plot") {
		plot(x, y, xlab = x.name, ylab = deparse(substitute(y)), \dots)
		lines(sm1$eval.points, sm1$estimate, col = 2)
	}
	sm1$aux <- list(h.first = h.first, first.sm = ym, acf = autocorr, 
		raw.residuals = r)
	invisible(sm1)
}}
}
\keyword{nonparametric}
\keyword{smooth}
\keyword{regression}
% Converted by Sd2Rd version 0.3-2.
