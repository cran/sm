\name{sm.pca}
\alias{sm.pca}

\title{
Smooth principal components analysis
}

\description{
This function calculates principal components in a manner which changes 
smoothly with a covariate.  The smooth eigenvalues and eigenvector
loadings can be plotted.  A permutation test of equality of the components,
both eigenvalues and eigenvectors, can be carried out.    
}

\usage{
sm.pca(x, Y, h, cor = TRUE, nperm = 100, pc = 1, ...)
}

\arguments{
\item{x}{
either a vector of covariate values or a list object which is the output of a previous call to \code{sm.pca}.  In the latter case, previously computed information is used to create plots and tests and the arguments \code{Y}, \code{h}, \code{cor} and \code{nperm} are not required.
}
\item{Y}{
a matrix of responses whose rows correspond to the covariate values.
}
\item{h}{
the smoothing parameter which controls the smoothness of estimation with respect to the covariate \code{x}.
}
\item{cor}{
a logical value indicating whether the correlation, rather than covariance, matrix should be used.
}
\item{nperm}{
the number of permutations used in the permutation test and graphical reference band.
}
\item{pc}{
an integer value indicating the component to be plotted against the covariate.
}
\item{\dots}{
  other optional parameters are passed to the \code{sm.options}
  function, through a mechanism which limits their effect only to this
  call of the function. Those relevant for this function are the following:
  \code{display} (here set to \code{"eigevalues"} or \code{"eigenvectors"})
  \code{ngrid},
  \code{xlab};
see the documentation of  \code{\link{sm.options}} for their description.
}
}

\value{
a list with the following components:
\describe{
\item{xgrid}{a vector of values on the covariate scale at which the smooth estimates are constructed.}
\item{evals}{a matrix whose columns give the smooth eigenvalues for each component at the covariate values.}
\item{evecs}{a three-dimensional array whose third dimension corresponds to the covariate values and whose second dimension indexes the smooth components.}
\item{mhat}{a matrix whose columns give the estimated smooth means for each dimension of \code{Y} at the covariate values.}
\item{var.explained}{a matrix whose rows give the proportions of variance explained by each component at each covariate value.}
\item{xlab}{the label attached to the x-axis.}
\item{h}{the smoothing parameter used.}
\item{x}{the covariate values, after removal of missing cases.}
\item{Y}{the matrix of response values, after removal of missing cases.}
\item{cor}{a logical indicator of whether the correlation, rather than covariance, matrix is used in the construction of the eigenvalues and eigenvectors.}
}
When a test or reference band is computed, the list has the additional components:
\describe{
\item{nperm}{the number of permutations used.}
\item{evals.perm}{the eigenvalues computed from the permuted data.}
\item{evecs.perm}{the eigenvectors computed from the permuted data.}
}
When display contains \code{"eigenvalues"} or \code{"eigenvectors"}, the list has the additional components:
\describe{
\item{p.values}{the p-value for a test of constant eigenvalue for the component identified by \code{pc}.}
\item{p.vectors}{the p-value for a test of constant eigenvectors for the component identified by \code{pc}.}
}
When display contains \code{"eigenvalues"}, the list has the additional component:
\describe{
\item{band}{a matrix whose two columns contain the boundaries of a reference band which indicates where the smooth eigenvalue curve should like if the hypothesis of no change in the eigenvalues with the covariate is correct.}
}
When display contains \code{"eigenvectors"}, the list has the additional components:
\describe{
\item{xgrid.plot}{a vector of values used for plotting the smooth eigenvectors.}
\item{evecs.plot}{a matrix whose rows contain the smooth eigenvectors at each value of \code{xgrid.plot}.}
\item{evecs.plot}{a matrix whose columns contain the colours for the line segments in each smooth eigenvector component.}
}
}

\section{Side Effects}{
a plot on the current graphical device is produced, unless \code{display="none"}
}

\details{
Several further arguments may be set and these are passed to \code{sm.options}.  Relevant arguments for this function are \code{display} (\code{"eigenvalues"}, \code{"eigenvectors"}), \code{ngrid} and \code{df}.  See \code{link{sm.options}} for further details.

The smoothing is performed by the local constant kernel method and the smoothing parameter corresponds to the standard deviation of a normal kernel function.  If \code{h} is left unspecified then it is selected to correspond to the degrees of freedom set by the parameter \code{df}.

The reference band for a constant eigenvalue is constructed from the upper and lower pointwise 2.5 percentiles of the smooth eigenvalue curves from the data with permuted covariate values.  The p-value compares the observed value of the difference between the smoothed and constant eigenvalues, summed over the covariate grid in \code{eval.points}, with the values generated from the permuted data.

In the eigenvector case, a reference band is computed from the percentiles of the curves from the permuted data, for each of the loadings.  In order to plot all the loadings curves simultaneously, the locations where each curve lies inside its respective reference band are indicated by pale colour.  The p-value compares the observed value of \code{1 - sum(e*e0)^2}, where \code{e} and \code{e0} are the eigenvectors under the smooth and constant scenarios (summed over the covariate grid), with the values generated from the permuted data.  This test statistic differs from the one described in the Miller and Bowman (2012) reference below.  It has been used as it conveniently handles the arbitrary sign of principal components.

When some components explain similar proportions of variance, the eigenvalues and eigenvectors can easily interchange, causing apparent sharp changes in the eigenvalue and eigenvector curves.  It is difficult to track the components to avoid this.
}

\references{
Miller, C. and Bowman, A.W. (2012). 
Smooth principal components for investigating changes in covariances over time.
\emph{Applied Statistics} \bold{61}, 693--714.
}

\seealso{
\code{\link{sm.regression}}, \code{\link{sm.options}}
}

\examples{
\dontrun{
Y    <- log(as.matrix(aircraft[ , -(1:2)]))
year <- aircraft$Yr
h    <- h.select(year, Y[ , 1], method = "df", df = 4)
spca <- sm.pca(year, Y, h, display = "none")
sm.pca(year, Y, h, display = "eigenvalues")
sm.pca(year, Y, h, display = "eigenvectors", ylim = c(-1, 1))

# The following code shows how the plots can be redrawn from the returned object

spca <- sm.pca(year, Y, h, display = "eigenvalues")
spca <- sm.pca(year, Y, h, display = "eigenvectors", ylim = c(-1, 1))

with(spca, {
   ylim <- range(evals[ , 1], band)
   plot(xgrid, evals[ , 1], type = "n", ylab = "Variance", ylim = ylim)
   polygon(c(xgrid, rev(xgrid)), c(band[ , 1], rev(band[ , 2])),
           col = "lightgreen", border = NA)
   lines(xgrid, evals[ , 1], col = "red")
})

with(spca, {
   pc <- 1
   plot(range(xgrid.plot), range(evecs.plot), type = "n",
        xlab = "x", ylab = "PC loadings")
   for (i in 1:ncol(Y))
      segments(xgrid.plot[-length(xgrid.plot)],
               evecs.plot[-nrow(evecs.plot), i],
               xgrid.plot[-1], evecs.plot[-1, i],
               col = col.plot[ , i], lty = i)
})
}
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
