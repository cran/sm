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
the value of the smoothing parameter located by the Sheather-Jones method.
}
\details{
See Section 2.4.4 of the reference below.
}

\note{As from version 2.1 of the package, a similar effect can be
  obtained with the new function \code{h.select}, via
  \code{h.select(x, method="sj")}. Users are encouraged to adopt
  this route, since \code{hsj} might be not accessible directly
  in future releases of the package.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). \emph{Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
   \code{\link{h.select}}, \code{\link{hnorm}}, \code{\link{hcv}}
}
\examples{
x <- rnorm(50)
hsj(x)
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
