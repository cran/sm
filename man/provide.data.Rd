\name{provide.data}
\alias{provide.data}
\title{
Making data available as data.frame
}
\description{
This function is no longer available in the sm package.  It should be replaced by the use of \code{attach}, if necessary.  Each dataset now also has its own help file.

It was a utility function, widely used in the scripts accompanying
the book described below.  The function provided access to the dataset identified by \code{name}.  For flexibility, the datasets were provided in ASCII form, with the name of each variable listed in the first row of the file.  This function reads the files and makes the data available as a data frame.
}
\usage{
provide.data(data, path, options = list())
}
\arguments{
\item{data}{
name of the data to be loaded and attached as \code{data.frame}
}
\item{path}{
the path where the data and its documentation should be searched for,
The default value is an appropriate sub-directory of the \code{sm} package. 
}
\item{options}{
  A list of options passed to \code{\link{sm.options}}.  The one used is
  \code{describe}, a logical flag.  If \code{describe=TRUE} (default), a
  documentation file of the data is searched and printed, if available.
}}
\value{
none
}
\section{Side Effects}{
messages are printed on the command window, describing  progress of
the operation. If \code{describe=TRUE} and a documentation file exists, this
is printed on the command windows or another windows, depending on
the type of platform where the program is executed.
}
\details{
the data file is assumed to be called \code{data}.dat and the documentation
file describing the data (if present) is assumed to be called \code{data}.doc.
If the \code{data.frame} is already attached, it is re-attached in the second
position of the \code{search} list. 

To set \code{describe=FALSE} for the rest of the current session,  use 
\code{sm.options(describe=FALSE)}

The function can easily be adapted to play a similar role for
other packages.
}
\author{
Bowman, A.W. and Azzalini, A.
}
\seealso{
  \code{\link{data.frame}}, \code{\link{attach}}, \code{\link{sm}},
  \code{\link{sm.options}}
}
\examples{
provide.data(birth)
}
\keyword{utilities}
% Converted by Sd2Rd version 1.15.
