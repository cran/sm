\name{sm.surface3d}
\alias{sm.surface3d}
\title{
Adding a regression surface to an rgl plot.
}
\description{
This function adds a regression surface, defined by a matrix of heights
at a regular grid of values of two covariates, to an \code{rgl} plot.
Missing values can be accommodated.
}
\usage{
 sm.surface3d(eval.points, surf, scaling, 
                    col = "green", col.mesh = "black", 
                    alpha = 0.7, alpha.mesh = 1, lit = TRUE, \dots)
}
\arguments{
\item{eval.points}{
if this is a two-column matrix then each column defines the marginal grids of covariate values.  Alternatively, a list with two components can also be used to handle cases where the grids are of different size.
}
\item{surf}{
a matrix of heights corresponding to the grid of covariate values.
NAs are allowed.
}
\item{scaling}{
a function to define the scaling for the \code{rgl} plot.  This
function is returned by an initial call to \code{rp.plot3d} in the
\code{rpanel} package.
}
\item{col}{
the colour of the surface.  If \code{col} is set to a single value, this
is replicated across the two components.  However, a matrix of values
corresponding to the entries of \code{surf} can also be supplied.
}
\item{col.mesh}{
the colour of the surface mesh.  If \code{col.mesh} is set to a single value, this
is replicated across the two components.  However, a matrix of values
corresponding to the entries of \code{surf} can also be supplied.
}
\item{alpha}{
the transparency of the filled triangles defining the surface.  Setting 
this to \code{0} will remove the filled triangles from the plot.
}
\item{alpha.mesh}{
the transparency of the lines drawn across the regular grid of covariate 
values.  Setting this to \code{0} will remove the lines from the plot.
}
\item{lit}{
a logical variable which controls whether the \code{rgl} plot is lit or not.
}
\item{\dots}{
other optional parameters which are passed to \code{material3d} in the 
\code{rgl} package.
}}
\value{
a vector of length 2 containing the ids of the filled surface and lines
added to the \code{rgl} plot.
}
\section{Side Effects}{
  a surface is added to the \code{rgl} plot.  
}
\details{
  the principal motivation for this function is that is can handle missing
  data in regression surfaces.  In particular, it can be used to plot the
  results of applying \code{sm.regression}.  In addition, the function can 
  be used to build up more complex plots by adding successive surfaces.
}
\seealso{
  \code{\link{sm.regression}}
}
\examples{
with(trawl, {
   Zone93    <- (Year == 1 & Zone == 1)
   Position  <- cbind(Longitude - 143, Latitude)
   model1 <- sm.regression(Position[Zone93,], Score1[Zone93],
        h= c(0.1, 0.1), display = "rgl", xlab="Longitude - 143")
   model2 <- sm.regression(Position[Zone93,], Score1[Zone93],
        h= c(0.2, 0.2), display = "none")
   sm.surface3d(model2$eval.points, model2$est, model1$scaling, col = "red")
})
}
\keyword{nonparametric}
\keyword{regression}
\keyword{smooth}
