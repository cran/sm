\name{provide.data}
\alias{provide.data}
\title{
Making data available as data.frame
}
\description{
This is a utility function, widely used in the scripts acconpanying
the book described below.
The function provides access to the dataset identified by \code{name}.
For flexibility, the datasets are provided in ASCII form, with the name
of each variable listed in the first row of the file.
This function reads the files and makes the data available as a data frame.
}
\usage{
provide.data(data, path, describe=T)
}
\arguments{
\item{data}{
name of the data to be loaded and attached as \code{data.frame}
}
\item{path}{
the path where the data and its documentation should be searched for,
The default value is an appropriate sub-directory of the \code{sm} package. 
}
\item{describe}{
logical flag; if \code{describe=T} (default), a documentation file of the data
is searched and printed, if available.
}}
\value{
none
}
\section{Side Effects}{
messages are printed on the command window, describing  progress of
the operation. If \code{describe=T} and a documentation file exists, this
is printed on the command windows or another windows, depending on
the type of platform where the program is executed.
}
\details{
the data file is assumed to be called \code{data}.dat and the documentation
file describing the data (if present) is assumed to be called \code{data}.doc.
If the \code{data.frame} is already attached, it is re-attached in the 2nd
position of the \code{search} list. 
}
\references{
Bowman, A.W. and Azzalini, A. (1997). Applied Smoothing Techniques for
Data Analysis: the Kernel Approach with S-Plus Illustrations.
Oxford University Press, Oxford.
}
\seealso{
\code{\link{data.frame}}, \code{\link{attach}}
}
\examples{
provide.data(birth)
}
\keyword{utilities}
% Converted by Sd2Rd version 0.3-2.
