\name{sig.trace}
\alias{sig.trace}
\title{
A significance trace for a hypothesis test
}
\description{
This function creates a significance trace for a hypothesis test based
on a nonparametric smoothing procedure.  The p-value of the test is
plotted against a range of smoothing parameters.
}
\usage{
sig.trace(expn, hvec, ...)
}
\arguments{
\item{expn}{
an S-Plus expression which should define the hypothesis test to be
performed, with the value of the smoothing parameter \code{h} omitted from
the function call.
}
\item{hvec}{
a vector of smoothing parameters for which the test will be performed.
}
\item{...}{
further arguments which will be passed to \code{sm.options}.
}}

\value{
a list containing vectors with the smoothing parameters and p-values.
}
\section{Side Effects}{
If the largest p-value is greater than 0.05 then a horizontal line at
0.05 will be superimposed on any plot, for reference.
}
\details{
see Section 5.2 of the reference below.


Only tests involving a univariate smoothing parameter may be used.
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
x <- runif(50, 0, 1)
y <- 5*x^2 + rnorm(50)
sig.trace(sm.regression(x, y, model = "linear", display="none"), 
        hvec = seq(0.05, 0.3, length = 10))
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
