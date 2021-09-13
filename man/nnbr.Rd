\name{nnbr}
\alias{nnbr}
\title{
nearest neighbour distances from data in one or two dimensions
}
\description{
This function calculates the \code{k} nearest neighbour distance from each
value in \code{x} to the remainder of the data.  In two dimensions, Euclidean
distance is used after standardising the data to have unit variance in
each component.
}
\usage{
nnbr(x, k)
}
\arguments{
\item{x}{
the vector, or two-column matrix, of data.
}
\item{k}{
the required order of nearest neighbour.
}}
\value{
the vector of nearest neighbour distances.
}
\details{
see Section 1.7.1 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). \emph{Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
none.
}
\examples{
x  <- rnorm(50)
hw <- nnbr(x, 10)
hw <- hw/exp(mean(log(hw)))
sm.density(x, h.weights=hw)
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
