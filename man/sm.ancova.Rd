\name{sm.ancova}
\alias{sm.ancova}
\title{
Nonparametric analysis of covariance
}
\description{
This function allows a set of nonparametric regression curves to be
compared, both graphically and formally in a hypothesis test.  A reference
model, used to define the null hypothesis, may be either equality or 
parallelism.
}
\usage{
sm.ancova(x, y, group, h, model="none", band=TRUE, test=TRUE,
          h.alpha=2 * diff(range(x))/length(x),
          weights = as.integer(rep(1, length(x))), \dots)
}
\arguments{
\item{x}{
a vector of covariate values.
}
\item{y}{
a vector of response values.
}
\item{group}{
a vector of group indicators.
}
\item{h}{
the smoothing parameter to be used in the construction of each of the
regression curves.
}
\item{model}{
a character variable which defines the reference model.  The values
\code{"none"}, \code{"equal"} and \code{"parallel"} are possible.
}
\item{band}{
a logical flag controlling the production of a reference band for the 
reference model.  A band will be produced only in the case of two groups.
}
\item{test}{
a logical flag controlling the production of a formal test, using the 
reference model as the null hypothesis.
}
\item{h.alpha}{
the value of the smoothing parameter used when estimating the vertical
separations of the curves under the parallelism model.
}
\item{weights}{
case weights.
}
\item{\dots}{
other optional parameters are passed to the \code{sm.options} function, through
a mechanism which limits their effect only to this call of the function; 
those relevant for this function are the following:
\itemize{
\item{display}{
any character setting other than \code{"none"} will cause a plot of the curves,
distinguished by line type, to be produced.
}
\item{ngrid}{
the size of the grid used to plot the curves.
}
\item{eval.points}{
a vector of points at which reference bands will be evaluated.
}
\item{xlab}{
the label attached to the x-axis.
}
\item{ylab}{
the label attached to the y-axis.
}}}}
\value{
a list containing an estimate of the error standard deviation and, where
appropriate, a p-value and reference model.  If the parallelism model has
been selected then a vector of estimates of the vertical separations of the 
underlying regression curves is also returned.
}
\section{Side Effects}{
a plot on the current graphical device is produced, unless \code{display="none"}
}
\details{
see Sections 6.4 and 6.5 of the book by Bowman \& Azzalini, and
the papers by Young \& Bowman listed below.
This function is a developed version of code originally written by Stuart Young.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.


Young, S.G. and Bowman, A.W. (1995).  
Nonparametric analysis of covariance.  
\emph{Biometrics}
\bold{51}, 920--931.


Bowman, A.W. and Young, S.G. (1996).  
Graphical comparison of nonparametric curves.  
\emph{Applied Statistics}
\bold{45}, 83--98.
}
\seealso{
\code{\link{sm.regression}}, \code{\link{sm.density.compare}}, \code{\link{sm.options}}
}
\examples{
x <- runif(50, 0, 1)
y <- 4*sin(6*x) + rnorm(50)
g <- rbinom(50, 1, 0.5)
sm.ancova(x, y, g, h = 0.15, model = "equal")
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
