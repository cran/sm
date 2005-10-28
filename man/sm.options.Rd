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
a list may be given as the only argument, or a vector of character
strings given as the only argument, or any number of arguments may be in
the \code{name=value} form , or no argument at all may be given.  See
the Value and Side Effects sections for explanation.

Arguments which are set by a function call will remain
in effect until the end of the
current S-plus session, unless overwritten by a subsequent call.  In
addition, they can be added as optional parameters of calls to specific
functions of the \code{sm} library; in this case, their effect is
limited to that function call.  See the documentation of specific
functions for the list of options which are recognised by that function.
Notice that some options are relevant only to some functions.
\itemize{
\item{add}{
a logical value which controls whether the density estimate is added to
the current plot.  Its default value is \code{FALSE}, which creates a new plot.
This argument applies only with one-dimensional data or to contours
created from two-dimensional data.
}
\item{band}{
 a logical value which controls whether the reference band appears on
 the plot. This argument applies only with one-dimensional data.
 Default: \code{band=FALSE}.
}
\item{col}{
the colour used for plotting observed points and estimated curves.
Where groups are used, \code{col} may be set to a vector of
colours associated with the groups.
Default: \code{col=1}.
}
\item{delta}{in \code{sm.density},
a value which will be added to the data before they are log transformed in
the procedure to handle positive data.  The value of \code{delta} is used
only when \code{positive} takes the value \code{TRUE}.  The default value
is the smallest value observed in each dimension.  This argument does not
apply with three-dimensional data.
Default: \code{delta=NA}
}
\item{describe}{
  logical flag whiich affects the behaviour of \code{sm.script} and
  \code{provide.data}. 
  If \code{describe=TRUE} (default), a data documentation file is printed. 
}
\item{df}{approximate degrees-of-freedom of the smoothing parameter used 
   in \code{sm.regression}, when a numerical value of \code{h} is not
   specified.
   In this case, the equivalent value of \code{h} will be computed and
   included in the list returned on exit from \code{sm.regression}.
   Default value is 6 if \code{x} is a vector and 12 if \code{x} is
   a matrix.
  }
\item{diff.ord}{in \code{sm.regression},
an integer defining the degree of differencing to be applied in the
estimation process.
When this argument is set to 1, the method of Rice,
based on the squared differences of pairs of neighbouring observations,
is used.  When the argument is set to 2 (deafult), the method of Gasser, 
Sroka and Jennen-Steinmetz, based on differences between each observation 
and a linear interplotation from its two neighbours, is used.
}
\item{display}{
This argument applies only with one- or two-dimensional data.  The
setting \code{"none"} will prevent any graphical output from being
produced.  In one dimensions, the default setting \code{"line"} will
produce the estimate, while the setting \code{"se"} will in
addition produce a variability band, showing the variability, but not
the bias, of the estimate.  In two dimensions, the default setting
\code{"persp"} will produce a perspective plot of the estimate,
while the settings \code{"slice"} and \code{"image"} will produce slice
(contour) or image plots.
}
\item{eval.grid}{
  logical flag which controls how the options \code{eval.points} are
  used fow two-dimensional data. If \code{eval.grid=TRUE} (default),
  evaluation is perfomed at points obtained by the cross-product of the
  two columns of \code{eval.points}. If \code{eval.grid=FALSE}
  then evaluation is performed at points with coordinates specified by
  the coordinates in \code{eval.points}
  }
\item{eval.points}{
the points at which the density or the regression curve or surface estimate 
should be evaluated, for the values returned in the result of the function.  
This should be a vector for one-dimensional data and a two-column matrix 
for two-dimensional data.
This argument does not apply with three-dimensional data.
}
\item{h.weights}{ 
a vector of weights which multiply the smoothing parameter used in the
kernel function at each observation.  This argument does not apply with 
three-dimensional data. Default value: 1.
}
\item{hmult}{
a factor which can be used to multiply the normal smoothing parameter
before construction of the density estimate. Default value: 1.
}
\item{hull}{  
a logical value which controls whether the estimate is evaluated and
plotted only on grid points which fall within the convex hull of the
data. When this argument is set to \code{FALSE}, evaluation and plotting
take place at all grid points where the contribution from at least one
kernel function is non-negligible.  Both of these settings ensure that
the estimate is not evaluated at points where there are no observations
nearby.  This argument applies only to \code{sm.regression} for the case
of two covariates.
}
\item{lty}{
  the line type used to plot the estimate.  This argument applies only
  when the estimate is displayed as a curve or a countour.
  Default value: 1.
}
\item{method}{
  the method used to select smoothing parameters.  In density estimation
  the default is \code{"normal"} which uses a vlue which is asymptotically
  optimal for the normal distribution.  Other possibilities are \code{"cv"}
  for cross-validation and \code{"sj"} for the Sheather-Jones method.
  In nonparametric regression, the
  deault is \code{"df"} which selects a smoothing parameters associated
  with the approximate degrees of freedom given in the \code{df} option.
  Other possibilities are \code{"cv"} for cross-validation and
  \code{"aicc"} for an AIC-based method proposed by Hurvich, Simonoff and
  Tsai.
}
\item{nbins}{
the number of bins used in one-dimensional binning operations;
in two-dimensional cases, \code{nbins} refers to the number of bins 
formed along each axis. Bins with 0 observed frequencies are ignored.
If \code{nbins=0}, binning is not performed; if \code{nbins=NA} (default),
binning is switched on when the number of observations exceeds
a certain threshold, which depends on the function.
}
\item{nboot}{
number of samples generated in bootstraps. Default value: 100.
}
\item{ngrid}{
the number of points in the regular grid used to plot the estimate.
For two- and three-dimensional data, \code{ngrid} refers to the
number of points along the axis in each dimension.
The same parameter is also used by a few other functions which perform some
form of search (e.g. \code{hcv}).
Default value for \code{sm.regression}:
  50 and 20 for 1-, 2-dimensional data, respectively.
Default value for \code{sm.density}:
  100, 50 and 20 for 1-, 2- and 3-dimensional data, respectively.
}
\item{panel}{
a logical value which controls the appearance of a panel display, in the
case of one- or two-dimensional data.  The panel can be used to alter the
value of the smoothing parameter.  For one-dimensional data it can also
launch an animation, and superimpose a normal reference band.  This argument 
does not apply with three-dimensional data.
}
\item{pch}{
  the standard plotting charater identified for data plotting.
  Default value: 1.
}
\item{phi}{
the vertical rotation (in degrees) of perspective plots of
estimate in the form of surfaces.  Default value: 40.
}
\item{poly.index}{
an integer defining local constant (0) or local linear (1) smoothing.
Default value: 1.
}
\item{positive}{
a logical value which indicates whether the data should be assumed to take
positive values only, in \code{sm.density}.
When this argument is set to \code{TRUE}, a log transformation
is applied to the data before construction of a density estimate.  The result
is transformed back to the original scale.  This argument does not apply with 
three-dimensional data. Default value: \code{FALSE}.
}
\item{props}{
a vector defining the proportions of the data to be included within each
contour in a slice plot, from two-dimensional data, or a contour surface
plot, from three-dimensional data.  In the three-dimensional case only
the first element of the vector will be used.  This argument does not apply 
to one-dimensional data.  Default value: \code{c(75,50,25)}.
}
\item{rugplot}{
logical flag which regulates whetehr a rugplot is superimposed to the
density estimate, in the univariate case. Dafault value: \code{TRUE}.
}
\item{show.script}{
  logical flag which affects the behaviour of \code{sm.script} when
  this is called with non-empty   argument. If \code{show.script=TRUE}
  (default) a window is opened to display the source code of the script.
  }
\item{structure.2d}{
  the structure of the smoothing parameter in two-dimensional settings.
  The default is \code{"scaled"}, which uses the structure
  (h*sd(x[,1]), h*sd(x[,2])).  Other possibilities are \code{"separate"},
  which uses (h1, h2), and \code{"common"} which uses (h, h).  The
  \code{"common"} option may be particularly appropriate when the data
  have a spatial origin, where distances in each variable havae the same
  meaning.  Note that the \code{"separate"} option is not available
  when \code{"method"} is set to \code{"df"}.
  }
\item{test}{
a logical flag controlling the production of a formal test, using the 
reference model as the null hypothesis. Default value: \code{TRUE}.
}
\item{theta}{
 the horizontal rotation (in degrees) of perspective plots of
 estimates in the form of surfaces.  Default value: -30.
}
\item{verbose}{
  regulates the amount of messages and other output printed out.
  If \code{verbose=0} only errors produce messages; if \code{verbose=1}
  (default value) warnings  and the more relevant numerical
  output are printed ; if \code{verbose=2} more messages and more
  numerical output are printed.
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
  the label attached to the z-axis (three-dimensional plots only).
}
\item{zlim}{ the range of the vertical axis when estimates are
  displayed as perspective plots.}
}}}
\value{
 a list with  the updated values of the parameters. If the argument list
 is not empty, the returned list is invisible.
}
\details{
  The parameter values set via a call to this function will remain in
  effect for the rest of the session, affecting the subsequent behaviour
  of the other functions of the \code{sm} library for which the given
  parameters are relevant.

  This offers the possibility of customizing the functioning of the \code{sm}
  library, for instance by insertions of appropriate calls to \code{sm.options}
  in a load hook for package \pkg{sm}.
}
\examples{
\dontrun{sm.options(poly.index = 0)
# subsequent regression estimations will be performed using local means
# instead of local regression
#
sm.options(describe = FALSE)  
# turns off typing documentation files of data loaded by `sm.script'
# (works from command-line)
# 
}}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
