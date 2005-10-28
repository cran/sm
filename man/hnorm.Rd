\name{hnorm}
\alias{hnorm}
\title{
Normal optimal choice of smoothing parameter in density estimation
}
\description{
This functions evaluates the smoothing parameter which is asymptotically
optimal for estimating a density function when the underlying distribution
is Normal.  Data in one, two or three dimensions can be handled.
}
\usage{
hnorm(x, weights)
}
\arguments{
\item{x}{
a vector, or matrix with two or three columns, containing the data.
}
\item{weights}{
  an optional vector of integer values
  which allows the kernel functions over the observations to take
different weights when they are averaged to produce a density estimate.  This
is useful, in particular, for censored data and to construct an estimate
from binned data. 
}}
\value{
the value of the asymptotically optimal smoothing parameter for Normal case.
}
\details{
See Section 2.4.2 of the reference below.
}
\note{As from version 2.1 of the package, a similar effect can be
  obtained with the new function \code{h.select}, via \code{h.select(x,
    method="normal", weights=weights)} or simply \code{h.select(x)}.
  Users are encouraged to adopt this route, since \code{hnorm} might be
  not accessible directly in future releases of the package.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
   \code{\link{h.select}}, \code{\link{hcv}}, \code{\link{hsj}}
}
\examples{
x <- rnorm(50)
hnorm(x)
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
