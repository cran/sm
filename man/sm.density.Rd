\name{sm.density}
\alias{sm.density}

\title{
Nonparametric density estimation in one, two or three dimensions.
}
\description{
This function creates a density estimate from data in one, two or three
dimensions.  In two dimensions a variety of
graphical displays can be selected, and in three dimensions a contour
surface can be plotted.  A number of other features of the construction
of the estimate, and of its display, can be controlled.

If the \code{rpanel} package is available, an interactive panel can be 
activated to control various features of the plot.

If the \code{rgl} package is also available, rotatable plots are
available for the two- and three-dimensional cases.  (For 
three-dimensional data, the \code{misc3d} package is also required.)

}
\usage{
sm.density(x, h, model = "none", weights = NA, group=NA, \dots)
}
\arguments{
\item{x}{a vector, or a matrix with two or three columns, containing the data.
         }
\item{h}{a vector of length one, two or three, defining the smoothing parameter.
         A normal kernel function is used and \code{h} is its standard deviation.
         If this parameter is omitted, a normal optimal smoothing parameter is used.
         }
\item{model}{This argument applies only with one-dimensional data.  Its default value
         is \code{"none"}.  If it is set to \code{"Normal"} (or indeed any value
         other than \code{"none"}) then a reference band, indicating where a
         density estimate is  likely to lie when the data are normally
         distributed, will be superimposed  on any plot.
         }
\item{weights}{a vector of integers representing frequencies of individual observations.
         Use of this parameter is incompatible with binning; hence \code{nbins} must 
         then be set to 0 or left at its default value \code{NA}.
         }
\item{group}{ a vector of groups indicators (numeric or character values) or
         a factor.
         }
\item{\dots}{
  other optional parameters are passed to the \code{sm.options} function,
  through a mechanism which limits their effect only to this call of the
  function. Those specifically relevant for this function are the following:
   \code{hmult}, 
   \code{h.weights}, 
   \code{band}, 
   \code{add}, 
   \code{lty}, 
   \code{display}, 
   \code{props}, 
   \code{xlab}, 
   \code{ylab}, 
   \code{zlab}, 
   \code{xlim}, 
   \code{ylim}, 
   \code{yht},        
   \code{nbins}, 
   \code{ngrid}, 
   \code{eval.points}, 
   \code{panel}, 
   \code{positive}, 
   \code{delta}, 
   \code{theta}, 
   \code{phi};
  see the documentation of  \code{\link{sm.options}} for their description.
}}
\value{
a list containing the values of the density estimate at the evaluation points,
the smoothing parameter, the smoothing parameter weights and the kernel 
weights.  For one- and two-dimensional data, the standard error of the estimate
(on the square root scale, where the standard error is approximately constant)
and the upper and lower ends of a variability band are also supplied.  Less 
information is supplied when the smoothing parameter weights
or kernel weights are not all 1, or when \code{positive} is set to \code{TRUE}.
}
\section{Side Effects}{
  a plot is produced, unless the option \code{display="none"} is set.
}
\details{
see Chapters 1, 2 and 6 of the reference below.
In the three-dimensional case, the contours of the density estimate are
constructed by the \code{contour3d} function in the \code{misc3d}
package of Feng & Tierney.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
  \code{\link{h.select}}, \code{\link{hnorm}}, \code{\link{hsj}}, \code{\link{hcv}},
  \code{\link{nise}}, \code{\link{nmise}}, \code{\link{sm}},
  \code{\link{sm.sphere}}, \code{\link{sm.regression}},
  \code{\link{sm.options}}
}
\examples{
#  A one-dimensional example
y <- rnorm(50)
sm.density(y, model = "Normal")
# sm.density(y, panel = TRUE)

#  A two-dimensional example
y <- cbind(rnorm(50), rnorm(50))
sm.density(y, display = "image")
# sm.density(y, panel = TRUE)


#  A three-dimensional example
# y <- cbind(rnorm(50), rnorm(50), rnorm(50))
# sm.density(y)
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
