\name{sm.sigma}
\alias{sm.sigma}
\title{
Estimation of the error standard deviation in nonparametric regression
}
\description{
This function uses ideas of local differencing to estimate the standard
deviation of the errors in a nonparametric regression model.  Simple
first-order differencing of pairs of neighbouring observations, or a method 
based on pseudo-residuals constructed from three neighbouring observations,
may be used.
}
\usage{
sm.sigma(x, y, diff.ord=2)
}
\arguments{
\item{x}{
a vector of covariate values.
}
\item{y}{
a vector of response values.
}
\item{diff.ord}{
an integer defining the degree of differencing to be applied in the
esimtation process.  When this argument is set to 1, the method of Rice,
based on the squared differences of pairs of neighbouring observations,
is used.  When the argument is set to 2, the method of Gasser, Sroka
and Jennen-Steinmetz, based on differences between each observation and
a linear interplotation from its two neighbours, is used.
}}
\value{
the estimated standard deviation of the errors.
}
\section{Side Effects}{
none.
}
\details{
see Section 4.3 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm.regression}}, \code{\link{sm.ancova}}
}
\examples{
x <- runif(100,-2, 2)
y <- x^2 + rnorm(50)
sm.sigma(x, y)
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 0.3-2.
