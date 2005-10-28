\name{sm.regression}
\alias{sm.regression}
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
 sm.regression(x, y, h, design.mat = NA, model = "none", weights =
                 NA, group = NA, \dots)
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
\code{"none"}, \code{"no effect"}, \code{"linear"} and \code{"parallel"}
are possible.
}
\item{weights}{
a vector  which allows the kernel functions associated with the observations 
to take different weights.  This is useful, in particular, when different
observations have different precisions.
The normal usage of this parameter is to associate observations with
frequencies; if the  \code{weights} are not integers, they are converted
to  integers, but in this case the standard errors and tests which are
computed cannot be considered.
This argument applies only to the case of one covariate.
Use of this parameter is incompatible with \code{binning}; hence
\code{nbins} must then be set to 0 or left at its default value \code{NA}.
}
\item{group}{
  a vector of groups indicators (numeric or character values) or
  a factor
}
\item{\dots}{
  other optional parameters are passed to the \code{sm.options}
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
    \code{xlab},
    \code{ylab},
    \code{zlab},
    \code{hull},
    \code{panel},
    \code{lty},
    \code{col};
  see the documentation of  \code{\link{sm.options}} for their description.
}}
\value{
a list containing the values of the estimate at the evaluation points,
the smoothing parameter and the smoothing parameter weights.
If a reference model has been specified and \code{test} set to
\code{TRUE}, then the p-value of the test
is also returned.  When there is only one covariate, the weights associated
with different obserations, an estimate of the error standard deviation and 
the standard error of the estimate are also returned.  If a reference model 
has been specified, this standard error refers to the comparison between 
the estimate and the reference model, and the values defining the reference 
model are also returned.
}
\section{Side Effects}{
  a plot on the current graphical device is produced, unless the option 
  \code{display="none"} is set.  
}
\details{
  see Chapters 3, 4 and 5 of the reference below.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis:}
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
  \code{\link{hcv}}, \code{\link{sm}}, \code{\link{sm.ancova}},
  \code{\link{sm.binomial}}, \code{\link{sm.poisson}},
  \code{\link{sm.regression.autocor}}, \code{\link{sm.survival}},
  \code{\link{sm.options}}
}
\examples{
#  An example with one covariate
x <- runif(100,-2, 2)
y <- x^2 + 2*x^3 + rnorm(50)/3
sm.regression(x, y, h=0.5)
sm.regression(x, y, col=2, add=TRUE)
#
u <-  sm.regression(x, y, h=h.select(x,y,method="cv"), display="none")
print(u$h)
#
#  An example with two covariates
x <- cbind(runif(100,-2, 2), runif(100,-2, 2))
y <- x[,1]^2 + x[,2]^2 + 0.7*x[,1]*x[,2] + rnorm(50)/3
sm.regression(x, y, h=c(1,1))
sm.regression(x, y, display="slice")
}
\keyword{nonparametric}
\keyword{regression}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
