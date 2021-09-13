\name{sm.regression}
\alias{sm.regression}
\title{Nonparametric regression with one or two covariates.}
\description{
This function creates a nonparametric regression estimate from data 
consisting of a single response variable and one or two covariates.
In two dimensions a perspective, image (\code{image}), contour (\code{slice}) 
or \code{rgl} plot of the estimated regression surface is produced.  
A number of other features of the construction of the estimate, and of 
its display, can be controlled.

If the \code{rpanel} package is available, an interactive panel can be activated
to control various features of the plot.  
}
\usage{
 sm.regression(x, y, h, design.mat = NA, model = "none", weights = NA,
                 group = NA, \dots)
}
\arguments{
   \item{x}{a vector, or two-column matrix, of covariate values.}
   \item{y}{a vector of responses.}
   \item{h}{a vector of length 1 or 2 giving the smoothing parameter.  A normal kernel
            function is used and \code{h} is its standard deviation.}
   \item{design.mat}{the design matrix used to produce \code{y} when these are assumed to be the 
            residuals from a linear model.}
   \item{model}{a character variable which defines a reference model.  The settings
            \code{"none"}, \code{"no effect"} and \code{"linear"} and all valid.  Note that
            when a \code{group} argument is used then \code{model} should be set to 
            \code{"equal"} or \code{"parallel"}, as befits an analysis of covariance model.}
   \item{weights}{a vector  which allows the kernel functions associated with the observations 
            to take different weights.  This is useful, in particular, when different
            observations have different precisions.
            The normal usage of this parameter is to associate observations with
            frequencies; if the  \code{weights} are not integers, they are converted
            to  integers, but in this case the standard errors and tests which are
            computed cannot be considered.
            This argument applies only to the case of one covariate.
            Use of this parameter is incompatible with \code{binning}; hence
            \code{nbins} must then be set to 0 or left at its default value \code{NA}.}
   \item{group}{a vector of groups indicators (numeric or character values) or
            a factor.  If this argument is used then the data are passed to the \code{sm.ancova}
            function.  See the details of the \code{model} argument above in that case.}
   \item{\dots}{other optional parameters are passed to the \code{sm.options}
            function, through a mechanism which limits their effect only to this
            call of the function; those relevant for this function are the following:
            \code{display},
            \code{hmult},
            \code{h.weights},
            \code{poly.index},
            \code{band},
            \code{add},
            \code{ngrid},
            \code{eval.points},
            \code{se},
            \code{se.breaks},
            \code{period},
            \code{xlab},
            \code{ylab},
            \code{zlab},
            \code{hull},
            \code{panel},
            \code{panel.plot},
            \code{lty},
            \code{col},
            \code{col.band},
            \code{col.mesh},
            \code{col.points},
            \code{col.palette};
            see the documentation of \code{\link{sm.options}} for their description.}
   }
\value{a list containing the values of the estimate at the evaluation points,
   the smoothing parameter and the smoothing parameter weights.
   If a reference model has been specified and \code{test} set to
   \code{TRUE}, then the p-value of the test
   is also returned.  When there is only one covariate, the weights associated
   with different observations, an estimate of the error standard deviation and 
   the standard error of the estimate are also returned.  If a reference model 
   has been specified, this standard error refers to the comparison between 
   the estimate and the reference model, and the values defining the reference 
   model are also returned.
   If an \code{rgl} display is used, then the indices of the surface and lines
   used to create the display are returned.
   }
\section{Side Effects}{
   a plot on the current graphical device is produced, unless the option 
   \code{display="none"} is set.  
   }
\details{
   When \code{display} is set to \code{"persp"} or \code{"rgl"}, a number of 
   graphical options are available.  By setting the \code{col} parameter to
   \code{"height"} or \code{"se"}, the surface will be painted by colours to
   reinforce the perception of height or indicate the relative sizes of the
   standard errors respectively.  When \code{model} is not \code{"none"}, 
   the colour coding refers to the number of standard errors which separate
   the smooth regression surface and the nominated model at each position.
   The parameter \code{"se.breaks"}, whose default value is \code{c(-3, -2, 3, 3)} 
   can then be used to set the colour ranges.  In this case, \code{col.palette}
   must be set to a list of colours whose length is one greater than the length
   of the cut-points in \code{"se.breaks"}.  If this is not the case, the
   default colour palette 
   \code{rev(rainbow(length(opt$se.breaks) + 1, start = 0/6, end = 4/6))}.
   
   If the argument \code{col} is not set then surface painting will be determined
   by the setting of \code{se}.  If neither is set then colour painting will be
   activated by default if \code{model != "none"}. (In this latter case, the
   argument \code{band}, retained from earlier versions for compatibility, will
   also be examined.)
   
   When \code{display} is set to \code{"rgl"}, some additional parameters
   can be used to control details of the plot.  Transparency can be set by
   \code{alpha}, which lies between \code{0} and \code{1}.  When \code{alpha}
   is set to a vector of length two, the first component refers to the surface
   panels and the second to the surface mesh.  Setting a component of \code{alpha}
   to \code{0} will remove the corresponding feature from the plot.  \code{col.mesh},
   whose valid values match those of \code{col}, controls the colour of the surface
   mesh.  The logical parameter \code{lit} has the same meaning as in the \code{rgl}
   package; see \code{material3d}.
   
   When \code{panel} is set to \code{"TRUE"}, an interactive control panel is
   created if the \code{rpanel} package is available.
   
   If a covariate is on a cyclical scale, this can be incorporated by setting
   the \code{period} argument to a vector (of length 1 or 2) whose components give 
   the values of the periods, or NA if the covariate is not periodic.
   
   See Chapters 3, 4 and 5 of the first reference below for the details of the
   construction of the estimate and its standard error.  The second reference
   gives further details and examples of surface painting.
   }
\references{
   Bowman, A.W. and Azzalini, A. (1997). 
   \emph{Applied Smoothing Techniques for Data Analysis:}
   \emph{the Kernel Approach with S-Plus Illustrations.}
   Oxford University Press, Oxford.
   
   Bowman, A.W. (2006).
   Comparing nonparametric surfaces.
   \emph{Statistical Modelling}, 6, 279-299.

   }
\seealso{
   \code{\link{hcv}}, \code{\link{sm}}, \code{\link{sm.ancova}},
   \code{\link{sm.binomial}}, \code{\link{sm.poisson}},
   \code{\link{sm.regression.autocor}}, \code{\link{sm.survival}},
   \code{\link{sm.options}}, \code{\link{sm.surface3d}}
   }
\examples{
with(trawl, {
   Zone92   <- (Year == 0 & Zone == 1)
   Position <- cbind(Longitude - 143, Latitude)
   dimnames(Position)[[2]][1] <- "Longitude - 143"

   par(mfrow = c(2, 2))
   sm.regression(Longitude, Score1, method = "aicc", col = "red",
       model = "linear")
   sm.regression(Position[Zone92, ], Score1[Zone92], display = "image", 
       theta = 120)
   sm.regression(Position[Zone92, ], Score1[Zone92], df = 12, col = "se",
       theta = 120)
   sm.regression(Position[Zone92, ], Score1[Zone92], df = 12, col = "se", 
       model = "linear", theta = 120)
   par(mfrow = c(1, 1))
})

# sm.regression(Position[Zone92, 2:1], Score1[Zone92], display = "rgl", df = 12)
# sm.regression(Position[Zone92, 2:1], Score1[Zone92], display = "rgl", df = 12,
#       alpha = c(0.9, 1), col = "se", model = "linear")

# sm.regression(Position[Zone92, 1], Score1[Zone92], panel = TRUE)
# sm.regression(Position[Zone92,  ], Score1[Zone92], panel = TRUE)
# sm.regression(Position[Zone92,  ], Score1[Zone92], panel = TRUE, display = "rgl")

}
\keyword{nonparametric}
\keyword{regression}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
