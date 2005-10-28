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
sm.regression.autocor(x = 1:n, y, h.first, minh, maxh, method = "direct", \dots)
}
\arguments{
\item{y}{
vector of the response values
}
\item{h.first}{
the smoothing parameter used for the initial smoothing stage.
}
\item{x}{
  vector of the covariate values; if unset, it is assumed to
  be \code{1:length(y)}.
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
possible values are \code{"no.cor"}, \code{"direct"} (default),
and \code{"indirect"}.
}
\item{\dots}{
  other optional parameters are passed to the \code{sm.options}
  function, through a mechanism which limits their effect only to this
  call of the function. Those relevant for this function are the following:
  \code{ngrid},
  \code{display};
see the documentation of  \code{\link{sm.options}} for their description.
}}
\value{
a list as returned from sm.regression called with the new value of
smoothing parameter, with an additional term \code{$aux} added which contains
the initial value \code{h.first}, the estimated curve using \code{h.first}, 
the autocorrelation function of the residuals from the initial fit, 
and the residuals.
}
\section{Side Effects}{
a new suggested value for \code{h} is printed; also, if the parameter \code{display}
is not equal to \code{"none"}, graphical output is produced on the current 
graphical device.
}
\details{
see Section 7.5 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.regression}}, \code{\link{sm.autoregression}}
}
\keyword{nonparametric}
\keyword{smooth}
\keyword{regression}
% Converted by Sd2Rd version 1.15.
