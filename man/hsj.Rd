\name{hsj}
\alias{hsj}
\title{
Sheather-Jones choice of smoothing parameter for density estimation
}
\description{
This function uses the Sheather-Jones plug-in method of selecting a
smoothing parameter which is suitable for constructing a density estimate
in the one-dimensional case.
}
\usage{
hsj(x)
}
\arguments{
\item{x}{
a vector of data.
}}
\value{
the value of the smoothing parameter localted by the Sheather-Jones method.
}
\section{Side Effects}{
none.
}
\details{
see Section 2.4.4 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.
Oxford University Press, Oxford.
}
\seealso{
\code{\link{hnorm}}, \code{\link{hcv}}
}
\examples{
x <- rnorm(50)
hsj(x)
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 0.3-2.
