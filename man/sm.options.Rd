\name{sm.options}
\alias{sm.options}
\title{
Set or return options of sm library
}
\description{
Provides a means to control the behaviour of the \code{sm} library
such as the colour of the plotted lines, the size of the grid in 2-D
estimation, the set of evaluations points, and many others.
}
\usage{
sm.options(\dots)
}
\arguments{
\item{\dots}{
a list may be given as the only argument, or a vector of
character strings given as the only argument, or any number of 
arguments may be in the \code{name=value} form , or no argument at all may be given.
See the Value and Side Effects sections for explanation.

Setting of some argument will remain in effect until the end of the current 
S-plus session, unless overwritten by a subsequent call. In addition, they 
can be added as optional parameters of call to specific functions of the \code{sm} 
library; in this case, their effect is limited  to that function call. 
See the documentation of specific functions for the list of options 
which are recognised by that function.
}
\item{add}{
a logical value which controls whether the density estimate is added to
the current plot.  Its default value is \code{FALSE}, which creates a new plot.
This argument applies only with one-dimensional data.
}
\item{band}{
a logical value which controls whether the reference band appears on the plot. 
This argument applies only with one-dimensional data.  Default: \code{band=FALSE}.
}
\item{col}{
The colour used for plotting observed points and estimated curves.
Default: \code{col=1}.
}
\item{delta}{
a value which will be added to the data before they are log transformed in
the procedure to handle positive data.  The value of \code{delta} is used only
when \code{positive} takes the value \code{TRUE}.  The default value is the smallest value
observed in each dimension.  This argument does not apply with 
three-dimensional data. Default: \code{delta=NA} 
}
\item{describe}{
This affectcs the behaviour of \code{sm.script} and \code{provide.data}. 
If \code{describe=TRUE} (default), a data documentation file is printed. 
}
\item{diff.ord}{
an integer defining the degree of differencing to be applied in the
esimtation process.  When this argument is set to 1, the method of Rice,
based on the squared differences of pairs of neighbouring observations,
is used.  When the argument is set to 2 (deafult), the method of Gasser, 
Sroka and Jennen-Steinmetz, based on differences between each observation 
and a linear interplotation from its two neighbours, is used.


}
\item{display}{
This argument applies only with one- or two-dimensional data.  
The setting \code{"none"} will prevent any graphical output from being produced.  
In one dimension, the default setting \code{"estimate"} will produce the density
estimate, while the setting \code{"se"} will in addition produce a variability
band, showing the variability, but not the bias, of the estimate.  In two
dimensions, the default setting \code{"persp"} will produce a perspective plot
of the density estimate, while the settings \code{"slice"} and \code{"image"} will
produce slice (contour) or image plots.
}
\item{eval.points}{
the points at which the density or the regression curve or surface estimate 
should be evaluated, for the values returned in the result of the function.  
This should be a vector for one-dimensional data and a two-column matrix 
for two-dimensional data. 
This argument does not apply with three-dimensional data.
}
\item{eye.mult}{
a vector of length 3 defining the eye position for the perspective plot
in the two covariate case.  The interpretation of this argument is similar
to the argument \code{eye} in the \code{persp} function.  Here the argument
is automatically scaled to adjust for the different axis scalings.
}
\item{h.weights}{
a vector of weights which multiply the smoothing parameter used in the
kernel function at each observation.  This argument does not apply with 
three-dimensional data. Default: \code{h.weights=1}.
}
\item{hmult}{
a factor which can be used to multiply the normal smoothing parameter
before construction of the density estimate. Default: \code{h.mult=1}.
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
\item{lty}{
the line type used to plot the estimate.  This argument applies only with 
one covariate.
}
\item{nbins}{
The number of bins used in one-dimensional binning operations;
in two-dimensional cases, \code{nbins} refers to the number of bins 
formed along each axis. Bins with 0 observed frequencies are ignored.
If \code{nbins=0}, binning is not performed; if \code{nbins=NA} (default),
binning is switched on when the number of observations exceeds
a certain threshold, which depends on the function.
}
\item{ngrid}{
The number of points in the regular grid used to plot the estimate.
For two- and three-dimensional data, \code{ngrid} refers to the number of points
along the axis in each dimension.
The same parameter is also used by a few other functions which perform some
form of search (e.g. \code{hcv}).
}
\item{panel}{
a logical value which controls the appearance of a panel display, in the
case of one- or two-dimensional data.  The panel can be used to alter the
value of the smoothing parameter.  For one-dimensional data it can also
launch an animation, and superimpose a normal reference band.  This argument 
does not apply with three-dimensional data.
}
\item{pch}{
the standard plottig charater identfier for data plotting. Default: \code{pch=1}.
}
\item{phi}{
the vertical rotation (in radians) of the contour plot of the density
estimate from three-dimensional data.  The default value is \eqn{\pi/4}{pi/4}.
}
\item{poly.index}{
an integer defining local constant (0) or local linear (1) smoothing.
(default value is 1)
}
\item{positive}{
a logical value which indicates whether the data should be assumed to take
positive values only.  When this argument is set to \code{TRUE}, a log transformation
is applied to the data before construction of a density estimate.  The result
is transformed back to the original scale.  This argument does not apply with 
three-dimensional data. Default: \code{positive=FALSE}.
}
\item{props}{
a vector defining the proportions of the data to be included within each
contour in a slice plot, from two-dimensional data, or a contour surface
plot, from three-dimensional data.  In the three-dimensional case only
the first element of the vector will be used.  This argument does not apply 
to one-dimensional data.  Default is \code{c(75,50,25)}.
}
\item{rugplot}{
Logical flag which regulates whetehr a rugplot is superimposed to the
density estimate, in the univariate case. Dafault: \code{rugplot=TRUE}.
}
\item{test}{
a logical flag controlling the production of a formal test, using the 
reference model as the null hypothesis. Defaul: \code{test=TRUE}.
}
\item{theta}{
the horizontal rotation (in radians) of the contour plot of the density
estimate from three-dimensional data.  The default value is \eqn{\pi/4}{pi/4}.
}
\item{xlab}{
the label attached to the x-axis.
}
\item{xlim}{
the range of the horizontal axis of the plot.  This argument does not apply
with three-dimensional data.
}
\item{yht}{
the upper limit of the vertical axis in a plot of a one-dimensional density
estimate.  The lower limit is always set to 0.  This argument does not apply 
with two- or three-dimensional data.
}
\item{ylab}{
the label attached to the y-axis. 
}
\item{ylim}{
the range of the vertical axis of the plot.  This argument does not apply 
with three-dimensional data.
}
\item{zlab}{
the label attached to the z-axis (three-dimensional plots only)
}}
\value{
a list with  the updated values of the parameters; if the argument list
is not empty, the returned list is invisible.
}
\details{
the parameter values set via a call to this function will remain in effect
for the rest of the session, affecting the subsequent behaviour of the other 
functions of the \code{sm} library for which the given parameters are relevant.

This offers the possibility of customizing the functioning of the \code{sm}
library, for instance by insertions of approriate calls to \code{sm.options}
in \code{.First}. However, assignements of \code{sm.options} are done at frame 0,
hence this must be done using a call as shown in the example below.
}
\examples{
\dontrun{sm.options(poly.index=0)
# subsequent regression estimations will be performed using local means
# instead of local regression
#
sm.options(describe=F)  
# turns off typing documentation files of data loaded by `sm.script'
# (works from command-line)
# 
assign(".sm.Options", sm.options(describe = FALSE), envir=.GlobalEnv)
# this line can be inserted in `.First' for permanent customisation
}}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
