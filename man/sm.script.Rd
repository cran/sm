\name{sm.script}
\alias{sm.script}
\title{
Running a script associated to the sm library
}
\description{
This is a utility function to run scripts, usually those associated 
with the book described below.
}
\usage{
sm.script(name, path)
}
\arguments{
\item{name}{
the name of the file containing the code; a \code{.q} suffix will be appended
to the name. If \code{name} is missing, a list of the scripts associated to
the \code{sm} library will be displayed.
}
\item{path}{
the name of the path where to look for the script. If \code{path} is missing,
it is assumed to be the appropriate location for the scripts of the \code{sm}
library.
}}
\details{
This utility allows the illustrations of the reference below to be reproduced
easily, since each of them has an associated script.  The display of the
script file itself is controlled by the setting of the logical variable
\code{show.script} in \code{sm.options}.  This is set to \code{TRUE} by default.
}
\references{
Bowman, A.W. and Azzalini, A. (1997). 
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.
}
\seealso{
\code{\link{sm}},
}
\examples{
sm.script()
sm.script(speed)
}
\keyword{utilities}
% Converted by Sd2Rd version 1.15.
