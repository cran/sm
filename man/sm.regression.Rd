\name{sm.regression}
\alias{sm.regression}
\alias{sm.regression.1d}
\alias{sm.regression.2d}
\alias{sm.regression.eval.1d}
\alias{sm.regression.eval.2d}
\alias{sm.regression.test}
\title{
Nonparametric regression with one or two covariates.
}
\description{
This function creates a nonparametric regression estimate from data 
consisting of a single response variable and one or two covariates.
With one covariate, an interactive panel can be activated
to control the smoothing parameter.  In two dimensions a perspective
plot of the estimated regression surface is produced.  A number of other 
features of the construction of the estimate, and of its display, can be 
controlled.
}
\usage{
sm.regression(x, y, h, design.mat = NA, model = "none", test = TRUE, 
              weights = rep(1, nobs), \dots)
}
\arguments{
\item{x}{
a vector, or two-column matrix, of covariate values.
}
\item{y}{
a vector of reponses.
}
\item{h}{
a vector of length 1 or 2 giving the smoothing parameter.  A normal kernel
function is used and \code{h} is its standard deviation.
}
\item{design.mat}{
the design matrix used to produce \code{y} when these are assumed to be the 
residuals from a linear model.  
}
\item{model}{
a character variable which defines a reference model.  The values
\code{"none"}, \code{"no effect"} and \code{"linear"} are possible.
}
\item{test}{
a logical flag controlling the production of a formal test, using the 
reference model as the null hypothesis.
}
\item{weights}{
a vector which allows the kernel functions associated with the observations 
to take different weights.  This is useful, in particular, when different
observations have different precisions.  This argument applies only to the
case of one covariate.
Use of this parameter is incompatible with binning; hence \code{nbins} must 
then be set to 0 or left at its default value \code{NA}.
}
\item{\dots}{
other optional parameters are passed to the \code{sm.options} function, through
a mechanism which limits their effect only to this call of the function; 
those relevant for this function are the following:
\itemize{
\item{display}{
The setting \code{"none"} will prevent any graphical output from being produced.
With one covariate, the default setting \code{"lines"} will produce the density
estimate, while the setting \code{"se"} will in addition produce a variability
band, showing the variability, but not the bias, of the estimate.  In two
dimensions, the default setting \code{"persp"} will produce a perspective plot
of the estimate.
}
\item{hmult}{
a factor which can be used to multiply the smoothing parameter
before construction of the nonparametric regression estimate.
}
\item{h.weights}{
a vector of weights which multiply the smoothing parameter used in the
kernel function at each observation. 
}
\item{poly.index}{
an integer defining local constant (0) or local linear (1) smoothing.
}
\item{band}{
a logical flag controlling the production of a reference band for the 
reference model.  A band will be produced only when there is a single
covariate and two groups.
}
\item{add}{
a logical value which controls whether the estimate is added to
the current plot.  Its default value is \code{F}, which creates a new plot.
This argument applies only with one covariate.
}
\item{ngrid}{
the number of points in the regular grid used to plot the estimate.
With two covariates, \code{ngrid} refers to the number of points
along the axis in each dimension.
Default:  50 and 20 for 1-, 2-dimensional data, respectively.
}
\item{eval.points}{
the points at which the estimate should be evaluated, for the
values returned in the result of the function.  This should be a vector
for one covariate and a two-column matrix for two covariates.
}
\item{xlab}{
the label attached to the x-axis.
}
\item{ylab}{
the label attached to the y-axis.
}
\item{zlab}{
the label attached to the z-axis (two covariate case only).
}
\item{hull}{
a logical value which controls whether the estimate is evaluated and plotted
only on grid points which fall within the convex hull of the data.  When
this argument is set to \code{FALSE}, evaluation and plotting take place at all grid 
points where the contribution from at least one kernel function is 
non-negligible.  Both of these settings ensure that the estimate is not 
evaluated at points where there are no observations nearby.  This argument 
applies only to the case of two covariates.
}
\item{panel}{
a logical value which controls the appearance of a panel display, in the
case of one covariate.  The panel can be used to alter the value of the 
smoothing parameter and in particular it can launch an animation.
}
\item{lty}{
the line type used to plot the estimate.  This argument applies only with 
one covariate.
}
\item{col}{
an integer defining the colour used to plot the estimate.  This argument 
applies only with one covariate.
}
\item{eye.mult}{
a vector of length 3 defining the eye position for the perspective plot
in the two covariate case.  The interpretation of this argument is similar
to the argument \code{eye} in the \code{persp} function.  Here the argument
is automatically scaled to adjust for the different axis scalings.
}}}}
\value{
a list containing the values of the estimate at the evaluation points,
the smoothing parameter and the smoothing parameter weights.  If a reference
model has been specified and \code{test} set to \code{T}, then the p-value of the test
is also returned.  When there is only one covariate, the weights associated
with different obserations, an estimate of the error standard deviation and 
the standard error of the estimate are also returned.  If a reference model 
has been specified, this standard error refers to the comparison between 
the estimate and the reference model, and the values defining the reference 
model are also returned.
}
\section{Side Effects}{
a plot on the current graphical device is produced, unless \code{display="none"}.
}
\details{
see Chapters 3, 4 and 5 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{hcv}}, \code{\link{sm}}, \code{\link{sm.ancova}}, \code{\link{sm.binomial}}, \code{\link{sm.poisson}}, \code{\link{sm.regression.autocor}},
\code{\link{sm.survival}}, \code{\link{sm.options}}
}
\examples{
#  An example with one covariate
x <- runif(100,-2, 2)
y <- x^2 + rnorm(50)
sm.regression(x, y, h=0.5)


#  An example with one covariate
x <- cbind(runif(100,-2, 2), runif(100,-2, 2))
y <- x[,1]^2 + x[,2]^2 + rnorm(50)
sm.regression(x, y, h=c(1,1))
}
\keyword{nonparametric}
\keyword{regression}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
