\name{sm.density}
\alias{sm.density}
\alias{sm.density.1d}
\alias{sm.density.2d}
\alias{sm.density.3d}
\alias{sm.density.eval.1d}
\alias{sm.density.eval.2d}
\alias{sm.density.positive.1d}
\alias{sm.density.positive.2d}
\alias{sm.density.positive.grid}

\title{
Nonparametric density estimation in one, two or three dimensions.
}
\description{
This function creates a density estimate from data in one, two or three
dimensions.  In one dimension, an interactive panel can be activated
to control the smoothing parameter, in two dimensions a variety of
graphical displays can be selected, and in three dimensions a contour
surface can be plotted.  A number of other features of the construction
of the estimate, and of its display, can be controlled.
}
\usage{
sm.density(x, h, model = "none", weights = rep(1,nobs), \dots)
}
\arguments{
\item{x}{
a vector, or a matrix with two or three columns, containing the data.
}
\item{h}{
a vector of length one, two or three, defining the smoothing parameter.
A normal kernel function is used and \code{h} is its standard deviation.
If this parameter is omitted, a normal optimal smoothing parameter is used.
}
\item{model}{
This argument applies only with one-dimensional data.  Its default value
is \code{"none"}.  If it is set to \code{"Normal"} (or indeed any value other than
\code{"none"}) then a reference band, indicating where a density estimate is 
likely to lie when the data are normally distributed, will be superimposed 
on any plot.
}
\item{weights}{
a vector which allows the kernel functions over the observations to take
different weights when they are averaged to produce a density estimate.
Use of this parameter is incompatible with binning; hence \code{nbins} must 
then be set to 0 or left at its default value \code{NA}.
This argument does not apply with three-dimensional data.
}
\item{\dots}{
other optional parameters are passed to the \code{sm.options} function, through
a mechanism which limits their effect only to this call of the function; 
those relevant for this function are the following:
\itemize{
\item{hmult}{
a factor which can be used to multiply the normal smoothing parameter
before construction of the density estimate.
}
\item{h.weights}{
a vector of weights which multiply the smoothing parameter used in the
kernel function at each observation.  This argument does not apply with 
three-dimensional data.
}
\item{band}{
a logical value which controls whether the reference band appears on the plot.
Its default value is \code{TRUE}.  
This argument applies only with one-dimensional data.
}
\item{add}{
a logical value which controls whether the density estimate is added to
the current plot.  Its default value is \code{FALSE}, which creates a new plot.
This argument applies only with one-dimensional data.
}
\item{lty}{
the line type used to plot the density estimate.  This argument applies only 
with one-dimensional data.
}
\item{display}{
This argument applies only with one- or two-dimensional data.  The setting
\code{"none"} will prevent any graphical output from being produced.  In one
dimensions, the default setting \code{"estimate"} will produce the density
estimate, while the setting \code{"se"} will in addition produce a variability
band, showing the variability, but not the bias, of the estimate.  In two
dimensions, the default setting \code{"persp"} will produce a perspective plot
of the density estimate, while the settings \code{"slice"} and \code{"image"} will
produce slice (contour) or image plots.
}
\item{props}{
a vector defining the proportions of the data to be included within each
contour in a slice plot, from two-dimensional data, or a contour surface
plot, from three-dimensional data.  In the three-dimensional case only
the first element of the vector will be used.  This argument does not apply 
to one-dimensional data.  
}
\item{xlab}{
the label attached to the x-axis.
}
\item{ylab}{
the label attached to the y-axis.
}
\item{zlab}{
the label attached to the z-axis (three-dimensional plots only)
}
\item{xlim}{
the range of the horizontal axis of the plot.  This argument does not apply
with three-dimensional data.
}
\item{ylim}{
the range of the vertical axis of the plot.  This argument does not apply 
with three-dimensional data.
}
\item{yht}{
the upper limit of the vertical axis in a plot of a one-dimensional density
estimate.  The lower limit is always set to 0.  This argument does not apply 
with two- or three-dimensional data.           
}
\item{nbins}{
The number of bins used in one-dimensional binning operations;
in two-dimensional cases, \code{nbins} refers to the number of bins 
formed along each axis. Bins with 0 observed frequencies are ignored.
If \code{nbins=0}, binning is not performed; if \code{nbins=NA} (default),
binning is switched on when the number of observations divided by
the number of dimensions exceeds 500.
}
\item{ngrid}{
the number of points in the regular grid used to plot the estimate.
For two- and three-dimensional data, \code{ngrid} refers to the number of points
along the axis in each dimension. 
Default: 100, 50 and 20 for 1-, 2- and 3-dimensional data, respectively.
}
\item{eval.points}{
the points at which the density estimate should be evaluated, for the
values returned in the result of the function.  This should be a vector
for one-dimensional data and a two-column matrix for two-dimensional data.
This argument does not apply with three-dimensional data.
}
\item{panel}{
a logical value which controls the appearance of a panel display, in the
case of one- or two-dimensional data.  The panel can be used to alter the
value of the smoothing parameter.  For one-dimensional data it can also
launch an animation, and superimpose a normal reference band.  This argument 
does not apply with three-dimensional data.
}
\item{positive}{
a logical value which indicates whether the data should be assumed to take
positive values only.  When this argument is set to \code{TRUE}, a log
transformation
is applied to the data before construction of a density estimate.  The result
is transformed back to the original scale.  This argument does not apply with 
three-dimensional data.
}
\item{delta}{
a value which will be added to the data before they are log transformed in
the procedure to handle positive data.  The value of \code{delta} is used only
when \code{positive} takes the value \code{TRUE}.  The default value is
the smallest value
observed in each dimension.  This argument does not apply with 
three-dimensional data.
}
\item{theta}{
the horizontal rotation (in radians) of the contour plot of the density
estimate from three-dimensional data.  The default value is \eqn{\pi/4}{pi/4}.
}
\item{phi}{
the vertical rotation (in radians) of the contour plot of the density
estimate from three-dimensional data.  The default value is \eqn{\pi/4}{pi/4}.
}}}}
\value{
a list containing the values of the density estimate at the evaluation points,
the standard error of the estimate, the upper and lower ends of a variability 
band, the smoothing parameter, the smoothing parameter weights and the kernel 
weights.  Less information is supplied when the smoothing parameter weights
or kernel weights are not all 1, or when \code{positive} is set to \code{TRUE}.
}
\section{Side Effects}{
a plot on the current graphical device is produced, unless \code{display="none"}.
}
\details{
see Chapters 1, 2 and 6 of the reference below.
The code for drawing contours of a three-dimensional density estimate has 
been amended from an original function written by Stuart Young.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{hnorm}}, \code{\link{hsj}}, \code{\link{hcv}}, \code{\link{nise}}, \code{\link{nmise}}, \code{\link{sm}}, \code{\link{sm.sphere}}, \code{\link{sm.regression}},
\code{\link{sm.options}}
}
\examples{
#  A one-dimensional example
y <- rnorm(50)
sm.density(y, model = "Normal")


#  A two-dimensional example
y <- cbind(rnorm(50), rnorm(50))
sm.density(y, display = "image")


#  A three-dimensional example
y <- cbind(rnorm(50), rnorm(50), rnorm(50))
sm.density(y)
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
