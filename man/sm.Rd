\name{sm}
\alias{sm}
\title{
The sm package: summary information
}
\description{
This package implements nonparametric smoothing methods described
in the book of Bowman & Azzalini (1997)
}
\section{Main Features}{
The functions in the package use kernel methods to construct
nonparametric estimates of density functions and regression curves in a
variety of settings, and to perform some inferential operations.

Specifically, density estimates can be constructed for 1-, 2- and
3-dimensional data. Nonparametric regression for continuous data can be
constructed with one or two covariates, and a variety of
tests can be carried out.  Several other data types can be
handled, including survival data, time series, count and binomial
data.
}
\section{Functions}{
The main functions are \code{sm.density} and \code{sm.regression}; other
functions intended for direct access by the user are:
\code{h.select}, \code{binning},
\code{sm.ancova}, \code{sm.autoregression}, \code{sm.binomial},
\code{sm.binomial.bootstrap}, \code{sm.poisson}, \code{sm.poisson.bootstrap},
\code{sm.options}, \code{sm.rm}, \code{sm.script}, \code{sm.sphere},
\code{sm.survival}, \code{sm.ts.pdf}, \code{sm.variogram}, \code{sm.pca}.  There are undocumented functions which are called by these.
}
\section{Scripts}{
  The function \code{sm.script} is used to run a set of examples (called
  scripts) presented in the book quoted below. These scripts are
  associated with the package but the package can be used independently of
  them.  The scripts are generally based on the functions of the package
  \code{sm}, but a few of them make used of the \code{gam} package.
}
\details{
Missing data are allowed; they are simply removed, together with
the associated variates from the same case, if any.
Datasets of arbitrary size can be handled by the current version of  
\code{sm.density},  \code{sm.regression} and \code{sm.ancova}, using
binning operations. 
}
\section{Requirements}{
R version >= 3.1.0. The \code{gam} package is used by
some of the scripts launched via \code{sm.script}, but it is not
used by the functions of this package.
}
\section{Version}{
This is version 2.2.
The most recent version of the package can be obtained from the CRAN archive.
}
\section{Details}{
  The book by Bowman and Azzalini (1997) provides more detailed and
  background information.  Algorithmic aspects of the software are
  discussed by Bowman & Azzalini (2003).  Differences between the first
  version of the package, described in the book, and the current one are
  summarized in the file \code{history.txt} which is distributed with
  the package.
}
\author{
Adrian Bowman (School of Mathematics and Statistics, University of Glasgow, UK) and
Adelchi Azzalini (Dept Statistical Sciences, University of Padua, Italy).
Please send comments, error reports, etc. to the authors.
}
\section{Acknowledgements}{
Important contributions to prototype versions of functions for some
specific techniques included here were made by a succession of students;
these include Stuart Young, Eileen Wright, Peter Foster, Angela Diblasi,
Mitchum Bock and Adrian Hines.  We are grateful for all these
interactions.  These preliminary version have been subsequently
re-written for inclusion in the public release of the package, with the
exception of the functions for three-dimensional density estimation,
written by Stuart Young.  We also thank Luca Scrucca who made useful
comments and who has ported the software to XLispStat.
We are particularly grateful to Brian Ripley for substantial help in
the production of installation files, the creation of MS-Windows
versions, initial porting of the software from S-Plus to R and for
maintaining the package on CRAN for several years.
}
\references{
Bowman, A.W. and Azzalini, A. (1997).
\emph{Applied Smoothing Techniques for Data Analysis: }
\emph{the Kernel Approach with S-Plus Illustrations.}
Oxford University Press, Oxford.

Bowman, A.W. and Azzalini, A. (2003).
Computational aspects of nonparametric smoothing
with illustrations from the \code{sm} library.
\emph{Computational Statistics and Data Analysis}, \bold{42}, 545--560.
}
\section{Licence}{
This package and its documentation are usable under the terms of the "GNU
General Public License", a copy of which is distributed with the package.
}
\keyword{nonparametric}
\keyword{smooth}
% Converted by Sd2Rd version 1.15.
