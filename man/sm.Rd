\name{sm}
\alias{sm}
\title{
The library sm: summary information
}
\description{
This library implements nonparametric smoothing methods described 
in the book of Bowman & Azzalini (1997)
}
\section{Main Features}{
The functions in the library use kernel methods to construct nonparametric 
estimates of density functions and regression curves in a variety of
settings, and to perform some inferential operations.

Specifically, density estimates can be performed for 1-, 2- and 3-dimensional  
data. Nonparametric regresion for continuous data can be constructed with
one or two covariates, and a variety of goodness-of-fit test for linear models
can be carried out. Many other data types can be handled; these include
survival data, time series, count and binomial data. 
}
\section{Functions}{
The main functions are \code{sm.density} and \code{sm.regression}; other
functions intended for direct access by the user are: \code{binning},
\code{sm.ancova}, \code{sm.autoregression}, \code{sm.binomial},
\code{sm.binomial.bootstrap}, \code{sm.poisson}, \code{sm.poisson.bootstrap},
\code{sm.options}, \code{sm.rm}, \code{sm.script}, \code{sm.sphere}, \code{sm.survival},
\code{sm.ts.pdf}.  There are undocumented functions which are called by the
above ones.
}
\details{
Missing data are allowed; they are simply removed, togeter with
the associated variates from the same case, if any.

Datasets of arbitrary size can be handled by the current version of  
\code{sm.density},  \code{sm.regression} and \code{sm.ancova}, using binning operations.
}
\section{REquirements}{
The library has been tested on S-plus 3.x, 4.0, 5.1
}
\section{Version}{
You are using version 2 (November 2000). 
The most recent version of the library can be obtained from either of 
the WWW pages:
\url{http://www.stats.gla.ac.uk/~adrian/sm}
\url{http://www.stat.unipd.it/~azzalini/Book\_sm}
}
\section{Manual}{
There is no manual except for on-line documentation.
The book by Bowman and Azzalini (1997) provides more detailed and 
background information.
Algorithmic aspects of the software are discussed by Bowman & Azzalini (2001).
Differences between the first version  of the library and the current one are 
summarized in the file \code{history.txt} which is distributed with the library.
}
\author{
Adrian Bowman (Dept Statistics, University of Glasgow, UK) and
Adelchi Azzalini (Dept Statistical Sciences, University of Padua, Italy).
Please send comments, error reports, etc. to the authors via the
abovementioned WWW page.
}
\section{Acknowledgements}{
Important contributions to prototype versions of functions for some specific 
techniques included here were made by a succession of students; these
include Stuart Young, Eileen Wright, Peter Foster, Angela Diblasi, 
Mitchum Bock and Adrian Hines. We are grateful for all these interactions.
These preliminary version have been subsequently re-written for inclusion
in the public release of the library, with the exception of the functions 
for three-dimensional density estimation, written by Stuart Young.
We also thank Luca Scrucca for useful remarks and Brian Ripley for substantial 
help in the production of  installation files, leading to much improved 
versions with respect to our original ones, and for tools to produce the 
MS-windows version starting from the Unix one.
}
\references{
Bowman, A.W. and Azzalini, A. (1997).
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.

Bowman, A.W. and Azzalini, A. (2001).
Computational aspects of nonparametric smoothing, 
with illustrations from the \code{sm} library.
To appear.
}
\section{Licence}{
This library and its documentation are usable under the terms of the "GNU 
General Public License", a copy of which is distributed with the package.
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
