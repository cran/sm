\name{geyser}
\alias{geyser}
\non_function{}
\usage{data(geyser)}
\title{Old Faithful Geyser Data}
\description{
A version of the eruptions data from the ``Old Faithful'' geyser in Yellowstone
National  Park,  Wyoming. This version comes from
Azzalini and Bowman (1990) and is of continuous measurement from August
1 to August 15, 1985.

Some nocturnal duration measurements were coded as 2, 3 or 4 minutes,
having originally been described as `short', `medium' or `long'.
}
\format{A data frame with 299 observations on 2 variables.
  \tabular{lll}{
    \code{duration}  \tab numeric  \tab Eruption time in mins \cr
    \code{waiting}   \tab numeric  \tab Waiting time to next eruption \cr
  }
}
\seealso{\code{\link{faithful}}}
\references{
  Azzalini, A. and Bowman, A. W. (1990) A look at some
  data on the Old Faithful geyser.  \emph{Applied Statistics}
  \bold{39}, 357--365.
}
\keyword{datasets}