\name{sm.poisson.bootstrap}
\alias{sm.poisson.bootstrap}
\title{
Bootstrap goodness-of-fit test for a Poisson regression model
}
\description{
This function is associated with \code{sm.poisson} for the underlying fitting
procedure.
It performs a Pseudo-Likelihood Ratio Test for the goodness-of-fit of
a standard parametric Poisson regression of specified \code{degree} in the
covariate \code{x}.
}
\usage{
sm.poisson.bootstrap(x, y, h, nboot=100, degree=1, \dots)
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
\item{nboot}{
number of bootstrap samples (default=100).
}
\item{degree}{
specifies the degree of the fitted polynomial in \code{x} on the logit scale
(default=1).
}
\item{...}{
additional parameters passed to \code{sm.poisson}
}}
\value{
a list containing the observed value of the Pseudo-Likelihood Ratio Test
statistic, its observed p-value as estimated via the bootstrap method,
and the vector of estimated dispersion parameters when this value is not 
forced to be 1.
}
\section{Side Effects}{
Graphical output representing the bootstrap samples is produced on 
the current graphical device. 
The estimated dispersion parameter, the value of the test statistic
and the observed significance level are printed.
}
\details{
see Section 5.4 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.poisson}}, \code{\link{sm.logit.bootstrap}}
}
\examples{
sm.poisson.bootstrap(exposure.time, N.events, 0.5, degree=2)
}
\keyword{nonparametric}
\keyword{smooth}
\keyword{htest}
% Converted by Sd2Rd version 0.3-2.
