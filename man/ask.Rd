\name{ask}
\alias{ask}
\title{
Ask data to be typed on the keyboard.
}
\description{
This function prints a message and prompts the user to type in data  
on the keyboard.
}
\usage{
ask(message="Type in datum")
}
\arguments{
\item{message}{
a character value, normally something describing the required input
}}
\value{
the object typed in.
}
\examples{
\dontrun{n <- ask("sample size")
nv <- ask("vector of sample sizes")
# in the second case, the user could reply eg with `c(5,10,30)'
}}
\keyword{utilities}
% Converted by Sd2Rd version 1.15.
