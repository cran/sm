\name{sm.script}
\alias{sm.script}
\title{Run or Lists Scripts in Package sm}
\usage{
sm.script(name, path)
}
\arguments{
 \item{name}{The name of a script file, without its extension. This can
     be given with or without quotes.}
 \item{path}{A path to the script.}
}
\description{
    With arguments, this runs the script \code{name} in directory
    \code{path}, or the \code{exec} subdirectory of the installed
    package directory. With no arguments it lists the scripts available.
}
\details{
    The scripts should have extension \code{.q}, and there should be a
    file \code{Index.doc} describing the scripts. On Windows the index
    file will be displayed in a \code{notepad} window unless
    \code{options("pager")} is set.
}
\value{
    None.
}
\example{
sm.script(rc_alter)
}
\keyword{utility}
