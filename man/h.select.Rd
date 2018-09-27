\name{h.select}
\alias{h.select}
\title{
Selection of the smoothing parameter
}
\description{
This function selects a smoothing parameter for density estimation
in one or two dimensions and for nonparametric regression with one
or two covariates.  Several methods of selection are available.
}
\usage{
 h.select(x, y = NA, weights = NA, group = NA, ...)
}
\arguments{
\item{x}{
  a vector, or two-column matrix.
  }
\item{y}{
a vector of reponses, in regression case.
}
\item{weights}{
  a vector of integers representing frequencies of individual
  observations.  Use of this parameter is incompatible with
  \code{binning}; hence \code{nbins} must then be set to \code{0}
  or left at its default value \code{NA}.
}
\item{group}{
  a vector of groups indicators (numeric or character values) or
  a factor
}
\item{\dots}{
  other optional parameters are passed to the \code{sm.options}
  function, through a mechanism which limits their effect only to this
  call of the function.  There are three which are relevant for this
  function, namely \code{method}, which specifies the method of
  smoothing parameter selection, \code{df}, which specifies the
  approximate degrees of freedom associated with the selected
  smoothing parameter, and \code{structure.2d} which determines
  the form of the smoothing parameters in the two-dimensional case.
  A full description of these arguments are given in the documentation
  of \code{\link{sm.options}}.
}}
\value{
  the value of the selected smoothing parameter.
}
\section{Side Effects}{
  none
}
\details{
  see the two references below for discussion of the methods of
  smoothing parameter selection.
  
  If the sample size is large, binning will be employed.  In the 
  case of \code{method = "cv"} the answer will therefore be
  different from that obtained through the function \code{hcv}
  where binning is not used.

  When the \code{group} argument is set, the chosen method of
  smoothing parameter selection is applied to each group and the
  value returned is the geometric mean of these.  This is intended
  for use in \code{\link{sm.density.compare}} and
  \code{\link{sm.ancova}}, where
  the same smoothing parameter is used for all groups so that
  the principal bias terms cancel when the estimates are compared.
}
\references{
  Bowman, A.W. and Azzalini, A. (1997). 
  \emph{Applied Smoothing Techniques for Data Analysis:}
  \emph{the Kernel Approach with S-Plus Illustrations.}
  Oxford University Press, Oxford.

  Hurvich, C.M., Simonoff, J.S. and Tsai, C.-L. (1998).
  Smoothing parameter selection in nonparametric regression
  using an improved Akaike information criterion.
  \emph{J. R. Statistic. Soc., Series B}, 60, 271-293.
  }
\seealso{ \code{\link{sm}},
  \code{\link{hcv}}, \code{\link{hsj}}, \code{\link{hnorm}}
}
\examples{
x <- rnorm(50)
h.select(x)
h.select(x, method = "sj")

x <- matrix(rnorm(100), ncol = 2)
h.select(x)
sm.density(x, method = "cv")

x <- rnorm(50)
y <- x^2 + rnorm(50)
h.select(x, y)
sm.regression(x, y, method = "aicc")

x <- matrix(rnorm(100), ncol = 2)
y <- x[,1]^2 + x[,2]^2 + rnorm(50)
h.select(x, y, method = "cv", structure.2d = "common")
sm.regression(x, y, df = 8)
}
\keyword{nonparametric}
\keyword{regression}
\keyword{smooth}
