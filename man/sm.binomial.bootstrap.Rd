\name{sm.binomial.bootstrap}
\alias{sm.binomial.bootstrap}
\title{
Bootstrap goodness-of-fit test for a logistic regression model. 
}
\description{
This function is associated with \code{sm.binomial} for the underlying fitting
procedure.  It performs a Pseudo-Likelihood Ratio Test for the
goodness-of-fit of a standard parametric logistic regression of specified
\code{degree} in the covariate \code{x}.
}
\usage{
sm.binomial.bootstrap(x, y, N = rep(1, length(x)), h, degree = 1,
        fixed.disp=FALSE, ...)
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
\item{N}{
a vector containing the binomial denominators.
If missing, it is assumed to contain all 1's.
}
\item{degree}{
specifies the degree of the fitted polynomial in \code{x} on the logit scale
(default=1).}
\item{fixed.disp}{if \code{TRUE}, the dispersion
  parameter is kept at value 1 across the simulated samples, otherwise
  the dispersion parameter estimated from the sample is used to generate
  samples with that dispersion parameter (default=\code{FALSE}).
}
\item{\dots}{
additional parameters passed to \code{\link{sm.binomial}}.
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
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.binomial}}, \code{\link{sm.poisson.bootstrap}}
}
\examples{
\dontrun{sm.binomial.bootstrap(concentration, dead, N, 0.5, nboot=50)}
}
\keyword{nonparametric}
\keyword{smooth}
\keyword{htest}
\keyword{models}
% Converted by Sd2Rd version 1.15.
